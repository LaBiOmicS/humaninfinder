"""Coding-evidence controls for a peptide window inside the 16S rRNA gene.

A reading frame located inside ``MT-RNR2`` can be conserved for either of two reasons: selection on
the peptide it would encode, or the constraint that keeps the host 16S rRNA functional. Sequence
conservation on its own does not distinguish them, because the whole gene is conserved. This module
implements the controls that do, over a cohort of genomes.

**This module does not use the MT-RNR2 profiles, the Pfam domain, or any other part of the locus
caller.** It aligns each genome's gene to the human rCRS reference bundled in
``data/conservation_reference.fasta`` with a pairwise aligner, and defines the peptide window by its
rCRS coordinates. That independence is what makes the result a control on the locus interpretation
rather than a restatement of it; ``tests/test_conservation.py`` asserts it.

Controls
--------
``pos3_share``
    Share of the variation that falls on the third codon position. A protein-coding frame under
    purifying selection concentrates variation there, where substitutions are mostly synonymous. A
    region constrained as RNA has no reason to favour any frame, so the expectation averaged over
    frames is ~1/3.
``aa_identity``
    Mean amino-acid identity to the reference peptide over aligned codons.

Each statistic is reported for the peptide frame, for the two shifted frames of the same window
(which encode nothing), and against a null built from every window of the same length across the
whole gene in all three frames. ``MT-CO1``, a real protein-coding gene in the same genomes under the
same alignment, is the positive control: the test must separate its true frame from its shifted ones.

Motif controls rank a sub-window (for example the ``CLLL`` motif) against every window of equal codon
length in the gene, and measure whether a segment downstream of an in-frame stop relaxes in the
genomes where it cannot be part of a translated peptide.

Caveats that belong with any output of this module: variation is counted per site without correcting
for phylogeny, which affects the window and its null alike; the null windows overlap, so percentiles
are descriptive and not a test; and a 21-codon window gives the codon-position test little power, so
a negative result means "no signal distinguishable from the rRNA background", not "no selection".
"""

from __future__ import annotations

import json
import statistics
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq

from .fasta_lint import open_fasta_text

DATA_DIR = Path(__file__).resolve().parent / "data"
DEFAULT_REFERENCE = DATA_DIR / "conservation_reference.fasta"

#: rCRS (NC_012920.1) 1-based inclusive coordinates of the bundled reference genes.
GENE_COORDS = {"MT-RNR2": (1671, 3229), "MT-CO1": (5904, 7445)}
#: rCRS 1-based inclusive window of the human Humanin reading frame, inside MT-RNR2.
HUMANIN_WINDOW = (2633, 2695)
#: Sub-windows reported by the motif control, as 1-based inclusive codon numbers in the peptide.
DEFAULT_MOTIFS = {"CLLL": (8, 11), "TSEIDLPVK": (13, 21)}

_MATCH, _MISMATCH, _OPEN, _EXTEND = 2, -3, -5, -2
_MIN_COVERAGE = 0.8

_ref: dict[str, str] = {}
_aligner: PairwiseAligner | None = None
_code_cache: dict[tuple[str, int], str] = {}


def load_reference(path: str | Path | None = None) -> dict[str, str]:
    """Reference gene sequences, keyed by gene name, from the bundled FASTA."""
    recs = SeqIO.parse(str(path or DEFAULT_REFERENCE), "fasta")
    return {r.id: str(r.seq).upper() for r in recs}


def _init(local: bool, ref_path: str | None = None) -> None:
    """Initialise the per-process reference and aligner (also the pool initialiser)."""
    global _ref, _aligner
    _ref = load_reference(ref_path)
    _aligner = PairwiseAligner()
    _aligner.mode = "local" if local else "global"
    _aligner.match_score, _aligner.mismatch_score = _MATCH, _MISMATCH
    _aligner.open_gap_score, _aligner.extend_gap_score = _OPEN, _EXTEND
    if not local:
        _aligner.end_gap_score = 0


def _aa(codon: str, table: int) -> str | None:
    """Translated residue of one codon, or None if it is not a plain ACGT triplet."""
    if len(codon) != 3 or any(c not in "ACGT" for c in codon):
        return None
    key = (codon, table)
    if key not in _code_cache:
        _code_cache[key] = str(Seq(codon).translate(table=table))
    return _code_cache[key]


def _project(ref: str, target: str) -> tuple[str, float]:
    """The base of ``target`` aligned to each position of ``ref`` ('-' where none is).

    A target with nothing to align -- empty, or made only of ``N`` and other non-ACGT symbols --
    projects as all gaps with score 0, and the coverage filter in ``analyse_cohort`` then drops it.
    Biopython raises on an empty sequence and returns an empty alignment list for a target that
    cannot be aligned at all, so both are handled here rather than surfacing as a traceback.
    """
    out = ["-"] * len(ref)
    if not target or not any(c in "ACGT" for c in target):
        return "".join(out), 0.0
    # Take the first optimal alignment from the iterator: len() would count every optimal
    # alignment, and on repetitive targets that count overflows a 64-bit int (OverflowError).
    aln = next(iter(_aligner.align(ref, target)), None)
    if aln is None:
        return "".join(out), 0.0
    for (r0, r1), (t0, t1) in zip(*aln.aligned, strict=True):
        out[r0:r1] = target[t0:t1]
    return "".join(out), aln.score


def project_genome(item: tuple[str, str, bool]) -> tuple[str, dict[str, str]]:
    """Project one genome onto every reference gene.

    Cohorts of different species are aligned locally against the whole genome on both strands.
    Same-species cohorts (``population=True``) are aligned globally to the rCRS-coordinate segment,
    which is both faster and safer when the sequences are nearly identical.
    """
    gid, seq, population = item
    seq = seq.upper()
    out = {}
    for gene, (a, b) in GENE_COORDS.items():
        ref = _ref[gene]
        if population:
            out[gene] = _project(ref, seq[max(0, a - 61) : b + 60])[0]
        else:
            both = (seq, str(Seq(seq).reverse_complement()))
            out[gene] = max((_project(ref, s) for s in both), key=lambda x: x[1])[0]
    return gid, out


@dataclass
class GeneData:
    """Per-site and per-codon summaries of one gene across a cohort.

    ``population`` switches the meaning of site variation: across species it is the fraction of
    genomes differing from the reference; within one species it is whether the site segregates at
    all, which is the quantity population studies report.
    """

    ref: str
    var: list[float | None]
    nt_frac: list[float | None]
    aa_id: list[float | None]
    aa_inv: list[bool]
    n: int

    @classmethod
    def build(cls, ref: str, projections: list[str], population: bool, table: int = 2) -> GeneData:
        from collections import Counter

        haplotypes = Counter(projections)
        length = len(ref)
        var: list[float | None] = [None] * length
        nt_frac: list[float | None] = [None] * length
        for p in range(length):
            if ref[p] not in "ACGT":
                continue
            total = differing = 0
            for hap, weight in haplotypes.items():
                if hap[p] in "ACGT":
                    total += weight
                    differing += weight * (hap[p] != ref[p])
            if total:
                frac = differing / total
                nt_frac[p] = frac
                var[p] = (1.0 if differing else 0.0) if population else frac
        aa_id: list[float | None] = [None] * length
        aa_inv: list[bool] = [False] * length
        for c in range(length - 2):
            ref_aa = _aa(ref[c : c + 3], table)
            if ref_aa is None:
                continue
            total = same = 0
            for hap, weight in haplotypes.items():
                hap_aa = _aa(hap[c : c + 3], table)
                if hap_aa is None:
                    continue
                total += weight
                same += weight * (hap_aa == ref_aa)
            if total:
                aa_id[c] = same / total
                aa_inv[c] = same == total
        return cls(ref=ref, var=var, nt_frac=nt_frac, aa_id=aa_id, aa_inv=aa_inv, n=len(projections))


def window_stats(gene: GeneData, start: int, frame: int, width: int) -> dict[str, Any]:
    """Statistics of the window ``[start, start + width)`` read in codon frame ``frame``."""
    pos_var = [0.0, 0.0, 0.0]
    for p in range(start, min(start + width, len(gene.ref))):
        if gene.var[p] is not None:
            pos_var[(p - frame) % 3] += gene.var[p]
    total = sum(pos_var)
    first = start + ((frame - start) % 3)
    codons = [c for c in range(first, start + width - 2, 3) if gene.aa_id[c] is not None]
    return {
        "pos_var": [round(x, 4) for x in pos_var],
        "pos3_share": pos_var[2] / total if total else None,
        "var_per_site": total / width,
        "aa_identity": statistics.mean(gene.aa_id[c] for c in codons) if codons else None,
        "aa_invariant_codons": sum(gene.aa_inv[c] for c in codons),
        "codons": len(codons),
    }


def scan_windows(gene: GeneData, width: int, step: int = 3) -> list[dict[str, Any]]:
    """Every window of ``width`` nt across the gene, in all three frames: the null distribution."""
    return [
        {"start": s, "frame": f, **window_stats(gene, s, f, width)}
        for s in range(0, len(gene.ref) - width + 1, step)
        for f in range(3)
    ]


def percentile_rank(value: float, null: list[float]) -> float:
    """Fraction of the null at or below ``value``."""
    return sum(v <= value for v in null) / len(null) if null else float("nan")


def coding_evidence(
    gene: GeneData,
    window: tuple[int, int],
    gene_start_rcrs: int,
    motifs: dict[str, tuple[int, int]] | None = None,
) -> dict[str, Any]:
    """The frame, shifted-frame and rRNA-null controls for one peptide window."""
    width = window[1] - window[0] + 1
    start = window[0] - gene_start_rcrs
    frame = start % 3
    by_frame = {f: window_stats(gene, start, f, width) for f in range(3)}
    rows = scan_windows(gene, width)
    null_p3 = [r["pos3_share"] for r in rows if r["pos3_share"] is not None]
    null_aa = [r["aa_identity"] for r in rows if r["aa_identity"] is not None]
    null_var = [r["var_per_site"] for r in rows if r["frame"] == 0]
    target = by_frame[frame]
    out: dict[str, Any] = {
        "window_rcrs": list(window),
        "frame": frame,
        "by_frame": {("peptide" if f == frame else f"shifted_{(f - frame) % 3}"): by_frame[f] for f in range(3)},
        "null_windows": len(null_p3),
        # ``is not None``: a share or identity of exactly 0 is a measurement and must be ranked
        "pos3_share_rank": (
            percentile_rank(target["pos3_share"], null_p3) if target["pos3_share"] is not None else None
        ),
        "null_pos3_share_median": statistics.median(null_p3) if null_p3 else None,
        "aa_identity_rank": (
            percentile_rank(target["aa_identity"], null_aa) if target["aa_identity"] is not None else None
        ),
        "null_aa_identity_median": statistics.median(null_aa) if null_aa else None,
        "var_per_site_rank": percentile_rank(target["var_per_site"], null_var),
        "codon_detail": [],
    }
    for k in range((width) // 3):
        c = start + 3 * k
        out["codon_detail"].append(
            {
                "codon": k + 1,
                "aa": _aa(gene.ref[c : c + 3], 2),
                "aa_identity": None if gene.aa_id[c] is None else round(gene.aa_id[c], 4),
                "nt_var_by_position": [
                    None if gene.nt_frac[c + j] is None else round(gene.nt_frac[c + j], 4) for j in range(3)
                ],
            }
        )
    if motifs:
        out["motifs"] = {
            name: motif_rank(gene, start + 3 * (lo - 1), frame, hi - lo + 1) for name, (lo, hi) in motifs.items()
        }
    return out


def motif_rank(gene: GeneData, start: int, frame: int, n_codons: int) -> dict[str, Any]:
    """Rank a sub-window against every window of the same codon length in the gene."""
    span = 3 * n_codons
    length = len(gene.ref)

    def summarise(s: int, f: int) -> tuple[float, float] | None:
        first = s + ((f - s) % 3)
        codons = [c for c in range(first, s + span - 2, 3) if gene.aa_id[c] is not None]
        if len(codons) < n_codons:
            return None
        return (
            statistics.mean(gene.aa_id[c] for c in codons),
            sum(gene.aa_inv[c] for c in codons) / len(codons),
        )

    null = [v for s in range(0, length - span + 1, 3) for f in range(3) if (v := summarise(s, f)) is not None]
    target = summarise(start, frame)
    if target is None or not null:
        return {"codons": n_codons, "aa_identity": None}
    ids, invs = [x[0] for x in null], [x[1] for x in null]
    return {
        "codons": n_codons,
        "aa_identity": target[0],
        "aa_identity_rank": percentile_rank(target[0], ids),
        "null_aa_identity_median": statistics.median(ids),
        "invariant_fraction": target[1],
        "invariant_fraction_rank": percentile_rank(target[1], invs),
        "fully_invariant_windows": sum(v == 1.0 for v in invs),
        "null_windows": len(null),
    }


def positive_control(gene: GeneData, window_width: int = 63) -> dict[str, Any]:
    """``MT-CO1`` in frame against out of frame: the test's response on a real coding gene."""
    rows = scan_windows(gene, window_width)
    in_frame = [r["pos3_share"] for r in rows if r["frame"] == 0 and r["pos3_share"] is not None]
    off_frame = [r["pos3_share"] for r in rows if r["frame"] != 0 and r["pos3_share"] is not None]
    in_aa = [r["aa_identity"] for r in rows if r["frame"] == 0 and r["aa_identity"] is not None]
    off_aa = [r["aa_identity"] for r in rows if r["frame"] != 0 and r["aa_identity"] is not None]
    return {
        "pos3_share_in_frame_median": statistics.median(in_frame) if in_frame else None,
        "pos3_share_off_frame_median": statistics.median(off_frame) if off_frame else None,
        "aa_identity_in_frame_median": statistics.median(in_aa) if in_aa else None,
        "aa_identity_off_frame_median": statistics.median(off_aa) if off_aa else None,
        "windows": len(rows),
    }


def downstream_relaxation(
    reference: str,
    projections: dict[str, str],
    window_start_in_gene: int,
    segment_codons: tuple[int, int],
    table: int = 2,
) -> dict[str, Any]:
    """Does a segment relax in the genomes where an upstream stop takes it out of any peptide?

    For every genome: amino-acid divergence of the segment, nucleotide divergence of the whole gene
    (a proxy for phylogenetic distance), and whether an in-frame stop precedes the segment. If the
    segment is conserved because it is peptide, it must drift once the peptide is truncated; if it is
    conserved because it is rRNA, it must not.

    The groups are not independent samples -- stop-bearing genomes cluster in a few lineages -- so
    the caller must report the number of independent origins alongside these numbers.
    """
    rows = []
    lo, hi = segment_codons
    for gid, proj in projections.items():
        if sum(c in "ACGT" for c in proj) < _MIN_COVERAGE * len(reference):
            continue
        nt_total = nt_diff = 0
        for p in range(len(reference)):
            if reference[p] in "ACGT" and proj[p] in "ACGT":
                nt_total += 1
                nt_diff += proj[p] != reference[p]
        if not nt_total:
            continue
        stop_before = any(
            _aa(proj[window_start_in_gene + 3 * k : window_start_in_gene + 3 * k + 3], table) == "*"
            for k in range(lo - 1)
        )
        aa_total = aa_diff = 0
        for k in range(lo - 1, hi):
            c = window_start_in_gene + 3 * k
            r, t = _aa(reference[c : c + 3], table), _aa(proj[c : c + 3], table)
            if r is None or t is None:
                continue
            aa_total += 1
            aa_diff += t != r
        if aa_total < hi - lo + 1:
            continue
        rows.append(
            {
                "id": gid,
                "nt_div_gene": nt_diff / nt_total,
                "aa_div_segment": aa_diff / aa_total,
                "stop_before_segment": stop_before,
            }
        )
    out: dict[str, Any] = {"segment_codons": list(segment_codons), "genomes": len(rows), "per_genome": rows}
    for key, want in (("with_stop", True), ("continuous", False)):
        group = [r for r in rows if r["stop_before_segment"] is want]
        ratios = [r["aa_div_segment"] / r["nt_div_gene"] for r in group if r["nt_div_gene"] > 0]
        out[key] = {
            "n": len(group),
            "mean_nt_div_gene": statistics.mean(r["nt_div_gene"] for r in group) if group else None,
            "mean_aa_div_segment": statistics.mean(r["aa_div_segment"] for r in group) if group else None,
            "median_ratio": statistics.median(ratios) if ratios else None,
            "identical_segment": sum(r["aa_div_segment"] == 0 for r in group),
        }
    return out


def analyse_cohort(
    fasta_path: str | Path,
    population: bool = False,
    reference_path: str | Path | None = None,
    window: tuple[int, int] = HUMANIN_WINDOW,
    motifs: dict[str, tuple[int, int]] | None = None,
    workers: int = 1,
    limit: int = 0,
) -> dict[str, Any]:
    """Run every control over a cohort of genomes and return one JSON-serialisable record."""
    if motifs is None:
        motifs = DEFAULT_MOTIFS
    items: list[tuple[str, str, bool]] = []
    with open_fasta_text(fasta_path) as handle:
        for i, rec in enumerate(SeqIO.parse(handle, "fasta")):
            if limit and i >= limit:
                break
            items.append((rec.id, str(rec.seq), population))
    if not items:
        raise ValueError(f"no sequences in {fasta_path}")

    ref_arg = str(reference_path) if reference_path else None
    _init(local=not population, ref_path=ref_arg)
    if workers > 1:
        with ProcessPoolExecutor(workers, initializer=_init, initargs=(not population, ref_arg)) as pool:
            projected = list(pool.map(project_genome, items, chunksize=32))
    else:
        projected = [project_genome(it) for it in items]

    out: dict[str, Any] = {
        "input": str(fasta_path),
        "genomes": len(items),
        "mode": "population" if population else "interspecies",
        "reference": str(reference_path or DEFAULT_REFERENCE),
        "genes": {},
    }
    host = "MT-RNR2"
    for gene_name, (gene_start, _end) in GENE_COORDS.items():
        ref = _ref[gene_name]
        min_cov = _MIN_COVERAGE * len(ref)
        kept = {gid: p[gene_name] for gid, p in projected if sum(c != "-" for c in p[gene_name]) >= min_cov}
        dropped = len(projected) - len(kept)
        if not kept:
            # Every genome fell below the coverage floor: the cohort has nothing to measure here.
            out["genes"][gene_name] = {
                "genomes_used": 0,
                "genomes_dropped_low_coverage": dropped,
                "min_coverage": _MIN_COVERAGE,
            }
            continue
        data = GeneData.build(ref, list(kept.values()), population=population)
        entry: dict[str, Any] = {"genomes_used": data.n}
        if dropped:
            entry["genomes_dropped_low_coverage"] = dropped
            entry["min_coverage"] = _MIN_COVERAGE
        if gene_name == host:
            entry["coding_evidence"] = coding_evidence(data, window, gene_start, motifs)
            seg = motifs.get("TSEIDLPVK")
            if seg:
                entry["downstream_relaxation"] = downstream_relaxation(ref, kept, window[0] - gene_start, seg)
        else:
            entry["positive_control"] = positive_control(data, window[1] - window[0] + 1)
        out["genes"][gene_name] = entry
    return out


def write_report(record: dict[str, Any], output_prefix: str) -> str:
    """Write the full record as JSON and return the path."""
    path = f"{output_prefix}_conservation.json"
    with open(path, "w") as f:
        json.dump(record, f, indent=1, sort_keys=True, default=str)
        f.write("\n")
    return path


def format_summary(record: dict[str, Any]) -> str:
    """Human-readable digest of a record, for the terminal."""
    lines: list[str] = []
    host = record["genes"].get("MT-RNR2", {})
    ce = host.get("coding_evidence")
    if not ce:
        dropped = host.get("genomes_dropped_low_coverage")
        if dropped:
            cov = host.get("min_coverage", _MIN_COVERAGE)
            return (
                f"no result: all {dropped} genome(s) covered less than {cov:.0%} of MT-RNR2.\n"
                "Nothing was measured. Check that the input holds whole mitochondrial genomes."
            )
        return "no host-gene result"
    lines.append(f"genomes used: {host['genomes_used']} ({record['mode']})")
    if host.get("genomes_dropped_low_coverage"):
        cov = host.get("min_coverage", _MIN_COVERAGE)
        lines.append(f"genomes dropped: {host['genomes_dropped_low_coverage']} (below {cov:.0%} coverage of MT-RNR2)")
    lines.append("")
    lines.append(f"{'frame':<12} {'pos3_share':>11} {'aa_identity':>12} {'invariant':>10}")
    for name, st in ce["by_frame"].items():
        p3 = "-" if st["pos3_share"] is None else f"{st['pos3_share']:.3f}"
        aa = "-" if st["aa_identity"] is None else f"{st['aa_identity']:.3f}"
        lines.append(f"{name:<12} {p3:>11} {aa:>12} {st['aa_invariant_codons']:>4}/{st['codons']:<5}")
    lines.append("")
    rank = ce.get("pos3_share_rank")
    med = ce.get("null_pos3_share_median")
    lines.append(
        f"pos3_share percentile within the 16S null: "
        f"{'-' if rank is None else f'{100 * rank:.0f}'} "
        f"(null median {'-' if med is None else f'{med:.3f}'}, {ce['null_windows']} windows)"
    )
    vr = ce.get("var_per_site_rank")
    lines.append(f"variation per site percentile: {'-' if vr is None else f'{100 * vr:.0f}'}")
    pc = record["genes"].get("MT-CO1", {}).get("positive_control")
    if pc and pc["pos3_share_in_frame_median"] is not None:
        lines.append(
            f"MT-CO1 positive control: pos3_share {pc['pos3_share_in_frame_median']:.3f} in frame "
            f"vs {pc['pos3_share_off_frame_median']:.3f} off frame"
        )
    for name, m in (ce.get("motifs") or {}).items():
        if m.get("aa_identity") is None:
            continue
        lines.append(
            f"motif {name}: aa identity {m['aa_identity']:.3f}, percentile "
            f"{100 * m['aa_identity_rank']:.0f} (null median {m['null_aa_identity_median']:.3f}); "
            f"{m['fully_invariant_windows']} fully invariant windows of equal length in the gene"
        )
    rel = host.get("downstream_relaxation")
    if rel and rel["with_stop"]["n"]:
        a, b = rel["with_stop"], rel["continuous"]
        lines.append("")
        lines.append(
            f"downstream segment (codons {rel['segment_codons'][0]}-{rel['segment_codons'][1]}): "
            f"aa divergence {a['mean_aa_div_segment']:.3f} in {a['n']} genomes with an upstream stop "
            f"vs {b['mean_aa_div_segment']:.3f} in {b['n']} continuous"
        )
        lines.append("  (stop-bearing genomes cluster in few lineages: report independent origins, not genome count)")
    return "\n".join(lines)
