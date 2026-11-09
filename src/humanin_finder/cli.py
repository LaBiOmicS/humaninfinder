"""Command line interface (CLI) for HumaninFinder."""

from __future__ import annotations

import csv
import hashlib
import json
import logging
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import click
from Bio import SeqIO

from . import __version__
from .conservation import analyse_cohort, format_summary, write_report
from .core import (
    MITOCHONDRIAL_TABLES,
    STATUS_COMPLETE,
    STATUS_INTERRUPTED,
    calculate_overlap,
    find_sorfs,
    get_gathering_threshold,
    guided_evolutionary_rescue,
    locate_humanin_loci,
    locus_gap_report,
    mdp_display_name,
    run_hmm_search_pyhmmer,
)
from .fasta_lint import LintOptions, lint_fasta, open_fasta_text
from .shlp import annotate_shlps

logger = logging.getLogger("humanin-finder")

LARGE_SEQUENCE_BP = 50_000
LOCUS_ZOOM_PAD_BP = 600

# Input policy for the pipeline: whole nuclear chromosomes and assemblies must pass the linter
PIPELINE_LINT_OPTIONS = LintOptions(
    alphabet="dna",
    max_input_bytes=64 * 1024**3,
    max_line_length=1_000_000_000,  # unwrapped (single-line) chromosome records
    max_sequence_length=1_000_000_000,
    max_total_sequence_length=50_000_000_000,
)

CSV_COLUMNS = [
    "sample_id",
    "name",
    "status",
    "seq",
    "length",
    "strand",
    "frame",
    "start_0based",
    "end_0based",
    "start_1based",
    "end_1based",
    "circular_wrap",
    "hmm_bitscore",
    "hmm_evalue",
    "hmm_supported",
    "synteny_score",
    "locus_profile",
    "locus_density",
    "start_offset_codons",
    "evidence_source",
    "evidence_tier",
    "reference_evidence_tier",
    "identity_to_reference",
    "genetic_code",
    "locus_mode",
    "table_1_pep",
]

STATUS_RANK = {STATUS_COMPLETE: 1.0, STATUS_INTERRUPTED: 0.5}


def _fmt(value: float | None, spec: str, missing: str = "NA") -> str:
    return missing if value is None else format(value, spec)


def _locus_score(r: dict[str, Any]) -> float | None:
    """Bits of the locus: MT-RNR2 profile score of a syntenic locus, else the PF15040 domain score."""
    return r.get("synteny_score") if r.get("synteny_score") is not None else r.get("hmm_bitscore")


def _file_sha256(path: Path) -> str | None:
    """SHA-256 of a file, read in chunks; None if it cannot be read."""
    try:
        h = hashlib.sha256()
        with open(path, "rb") as f:
            for chunk in iter(lambda: f.read(1 << 20), b""):
                h.update(chunk)
        return h.hexdigest()
    except OSError:
        return None


def write_run_provenance(
    output_prefix: str,
    fasta_path: Path,
    options: dict[str, Any],
    result_count: int,
    gap_rows: int = 0,
) -> dict[str, Any]:
    """Write ``{prefix}_run.json`` recording exactly how this run was invoked.

    Written so that a result can be traced back to the command that produced it: the published
    tables of a cohort are reproducible only if the subcommand and options are known (``scan`` and
    ``predict`` legitimately give different row counts on the same input).

    Timestamps are the only non-deterministic field; the rest is a function of the invocation.
    """
    record = {
        "tool": "humanin-finder",
        "version": __version__,
        "command": " ".join(sys.argv),
        "argv": list(sys.argv),
        "finished_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "input": {
            "path": str(fasta_path),
            "sha256": _file_sha256(Path(fasta_path)),
        },
        "options": options,
        "result_rows": result_count,
        # peptides whose start codon fell in a gap, so no frame could be read (see {prefix}_gaps.csv)
        "start_codon_gap_rows": gap_rows,
    }
    path = f"{output_prefix}_run.json"
    with open(path, "w") as f:
        json.dump(record, f, indent=1, sort_keys=True)
        f.write("\n")
    return record


def _csv_cell(value: Any, as_float: bool) -> str:
    """One CSV cell (see ``write_csv`` for the formatting rules)."""
    if value is None:
        return ""
    if isinstance(value, bool):
        return "True" if value else "False"
    if as_float and isinstance(value, (int, float)):
        return str(float(value))
    return str(value)


def write_csv(rows: list[dict[str, Any]], columns: list[str], path: str) -> None:
    """Write ``rows`` as CSV over ``columns`` with the standard library.

    This replaced ``pandas.DataFrame(rows).to_csv(index=False)`` and keeps its output, with one
    deliberate exception fixed in 1.0.0b1: pandas turned an integer column holding any ``None`` into
    float64, so ``start_offset_codons`` was written as ``-1.0``. A codon offset is an integer, so
    integer columns are now written as integers whether or not they have empty cells. Columns that
    hold real floats keep ``str(float)`` for every value. Otherwise: empty field for ``None``,
    ``True``/``False`` for booleans, ``str`` otherwise, minimal quoting, ``\n`` line endings.
    """
    # Per column: format as float only when the column holds a real float (bools are not numeric).
    as_float: dict[str, bool] = {}
    for col in columns:
        values = [r.get(col) for r in rows]
        present = [v for v in values if v is not None]
        numeric = bool(present) and all(isinstance(v, (int, float)) and not isinstance(v, bool) for v in present)
        has_float = any(isinstance(v, float) for v in present)
        as_float[col] = numeric and has_float
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(columns)
        for row in rows:
            writer.writerow([_csv_cell(row.get(col), as_float[col]) for col in columns])


def write_deterministic_outputs(results: list[dict[str, Any]], output_prefix: str) -> None:
    """Generate deterministic CSV, FASTA, GFF3, and BED tracks with FAIR attributes."""
    out_dir = Path(output_prefix).parent
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1. CSV
    csv_file = f"{output_prefix}_results.csv"
    write_csv(results, CSV_COLUMNS, csv_file)

    # 2. FASTA
    fasta_file = f"{output_prefix}_results.fasta"
    with open(fasta_file, "w") as f:
        for r in results:
            header = (
                f">{r['sample_id']}|{r.get('name', 'Humanin')}|{r['status']}|"
                f"coords_{r['start_1based']}-{r['end_1based']}|evidence_{r.get('evidence_source') or 'NA'}"
            )
            f.write(f"{header}\n{r['seq']}\n")

    # 3. GFF3 (features crossing the origin of a circular landmark end beyond its length)
    gff3_file = f"{output_prefix}_results.gff3"
    landmarks: dict[str, dict[str, Any]] = {}
    for r in results:
        if r.get("seq_length"):
            landmarks.setdefault(r["sample_id"], {"length": r["seq_length"], "circular": r.get("seq_circular", False)})
    used_ids: dict[str, int] = {}
    with open(gff3_file, "w") as f:
        f.write("##gff-version 3\n")
        f.write("#coordinate-system 1-based inclusive [start, end]\n")
        for seqid, meta in landmarks.items():
            f.write(f"##sequence-region {seqid} 1 {meta['length']}\n")
        for seqid, meta in landmarks.items():
            if meta["circular"]:
                f.write(f"{seqid}\tHumaninFinder\tregion\t1\t{meta['length']}\t.\t+\t.\tID={seqid};Is_circular=true\n")
        for r in results:
            strand_char = "+" if r.get("strand", 1) == 1 else "-"
            start_1 = int(r["start_1based"])
            end_1 = int(r["end_1based"])
            score = _locus_score(r)
            base_id = f"{r['sample_id']}_{r.get('name', 'Humanin')}_{start_1}_{end_1}_{strand_char}"
            used_ids[base_id] = used_ids.get(base_id, 0) + 1
            feature_id = base_id if used_ids[base_id] == 1 else f"{base_id}.{used_ids[base_id]}"
            attr = (
                f"ID={feature_id};Name={r.get('name', 'Humanin-like')};status={r['status']};"
                f"hmm_score={_fmt(r.get('hmm_bitscore'), '.1f')};"
                f"hmm_supported={str(bool(r.get('hmm_supported', False))).lower()};"
                f"synteny_score={_fmt(r.get('synteny_score'), '.0f')};"
                f"locus_density={_fmt(r.get('locus_density'), '.3f')};"
                f"evidence_source={r.get('evidence_source', 'NA')};"
                f"evidence_tier={r.get('evidence_tier', 'NA')};"
                f"reference_evidence_tier={r.get('reference_evidence_tier') or 'NA'};"
                f"translation_table={r.get('genetic_code', 2)};"
                f"model_version={__version__}"
            )
            # SO:0000236 ORF: coordinates span the reading frame (including the stop codon when complete)
            f.write(
                f"{r['sample_id']}\tHumaninFinder\tORF\t{start_1}\t{end_1}\t{_fmt(score, '.1f', '.')}"
                f"\t{strand_char}\t.\t{attr}\n"
            )

    # 4. BED (origin-crossing features are split into two blocks with the same name)
    bed_file = f"{output_prefix}_results.bed"
    with open(bed_file, "w") as f:
        f.write("#track name=HumaninFinder description='0-based half-open genomic intervals'\n")
        f.write("#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\n")
        for r in results:
            strand_char = "+" if r.get("strand", 1) == 1 else "-"
            # Fully named loci 1000; partial evidence 500 (position without peptide, or domain
            # without position)
            base_name = str(r.get("name", "")).split("_T")[0]
            score_int = 500 if base_name.endswith(("_position", "_domain")) else 1000
            name = f"{r.get('name', 'Humanin')}_{r['status']}"
            start, end = int(r["start_0based"]), int(r["end_0based"])
            seq_length = r.get("seq_length")
            blocks = [(start, end)]
            if seq_length and end > seq_length:
                blocks = [(start, seq_length), (0, end - seq_length)]
            for b_start, b_end in blocks:
                f.write(f"{r['sample_id']}\t{b_start}\t{b_end}\t{name}\t{score_int}\t{strand_char}\n")


def filter_non_redundant(
    candidates: list[dict[str, Any]],
    overlap_threshold: float = 0.5,
) -> list[dict[str, Any]]:
    """Greedy suppression of overlapping candidates, keeping syntenic, then HMM-supported, then highest-HMM-score ones.

    Candidates only suppress each other on the same strand and under the same genetic code:
    antisense ORFs and alternative translation tables are distinct hypotheses.
    """
    if not candidates:
        return []

    sorted_cands = sorted(
        candidates,
        key=lambda c: (
            1.0 if c.get("synteny_score") is not None else 0.0,
            1.0 if c.get("hmm_supported") else 0.0,
            float(c.get("hmm_bitscore") or 0.0),
        ),
        reverse=True,
    )

    accepted = []
    for cand in sorted_cands:
        is_redundant = False
        c_len = max(1, cand["end_0based"] - cand["start_0based"])
        for acc in accepted:
            if cand.get("genetic_code") != acc.get("genetic_code") or cand.get("strand") != acc.get("strand"):
                continue
            overlap = calculate_overlap(
                cand["start_0based"], cand["end_0based"], acc["start_0based"], acc["end_0based"]
            )
            if overlap / c_len > overlap_threshold:
                is_redundant = True
                break
        if not is_redundant:
            accepted.append(cand)

    return accepted


def parse_genetic_codes(table_arg: str | int | list[int]) -> list[int]:
    """Parse comma-separated genetic code tables into a validated list of unique integers.

    ``"auto"`` returns an empty list: the table is then chosen per record (``auto_genetic_code``).
    """
    if isinstance(table_arg, str) and table_arg.strip().lower() == "auto":
        return []
    if isinstance(table_arg, int):
        return [table_arg]
    if isinstance(table_arg, list):
        return [int(t) for t in table_arg]
    tables: list[int] = []
    for part in str(table_arg).split(","):
        part = part.strip()
        if not part:
            continue
        try:
            val = int(part)
            if not (1 <= val <= 33):
                raise ValueError
            if val not in tables:
                tables.append(val)
        except ValueError as err:
            raise click.BadParameter(f"Genetic code table must be an integer between 1 and 33, got: {part!r}") from err
    return tables or [2]


def auto_genetic_code(description: str, seq_len: int) -> int:
    """Vertebrate default code of a record: Table 2 for mitochondrial records, Table 1 otherwise.

    A record is mitochondrial when its header says so or when it is no longer than
    ``LARGE_SEQUENCE_BP`` (complete mitogenomes and mitochondrial contigs). Non-vertebrate
    mitogenomes use other codes (e.g. 5 for insects, 9 for echinoderms) and need ``--table``.
    """
    if "mitochondri" in description.lower() or seq_len <= LARGE_SEQUENCE_BP:
        return 2
    return 1


def _merge_regions(regions: list[tuple[int, int]]) -> list[tuple[int, int]]:
    merged: list[tuple[int, int]] = []
    for start, end in sorted(regions):
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def _make_candidate(
    rec_id: str,
    item: dict[str, Any],
    t_id: int,
    name_suffix: str,
    locus_mode: str,
) -> dict[str, Any]:
    """Build a result row from a sORF (find_sorfs) or a rescued peptide (guided_evolutionary_rescue)."""
    rescued = item.get("rescued", False)
    return {
        "sample_id": rec_id,
        "name_suffix": name_suffix,
        "status": item["status"],
        "seq": item["seq"],
        "dna_seq": item.get("dna_seq", ""),
        "length": item["length"],
        "strand": item["strand"],
        "frame": item["frame"],
        "start_0based": item["start"],
        "end_0based": item["end"],
        "start_1based": item["start"] + 1,
        "end_1based": item["end"],
        "circular_wrap": bool(item.get("circular_wrap", False)),
        "genetic_code": t_id,
        "locus_mode": locus_mode,
        "hmm_bitscore": item["hmm_score"] if rescued else None,
        "hmm_evalue": item["hmm_evalue"] if rescued else None,
        "synteny_score": item.get("synteny_score"),
        "locus_profile": item.get("locus_profile"),
        "locus_density": item.get("locus_density"),
        "vertebrate_like": item.get("vertebrate_like"),
        "start_offset_codons": item.get("start_offset_codons"),
        "identity_to_reference": item.get("identity_to_reference"),
        "table_1_pep": item.get("table_1_pep", ""),
    }


def _annotate_evidence(cand: dict[str, Any], gathering: float) -> None:
    """Derive the explicit evidence fields of a Humanin candidate (in silico homology only).

    A frame is reported when it lies at a Humanin locus located by its MT-RNR2 position or carries a
    PF15040 domain at the gathering threshold. The name says which routes hold: ``Humanin`` at a
    vertebrate-like locus or a locus frame with the domain, ``Humanin_position`` at a locus outside
    the vertebrate pattern, and ``Humanin_domain`` for a domain-bearing frame outside any locus
    (``core.mdp_display_name``). Frames without either route are not reported.
    """
    syn_ok = cand.get("synteny_score") is not None
    hmm_ok = cand.get("hmm_bitscore") is not None and cand["hmm_bitscore"] >= gathering
    base = mdp_display_name("Humanin", cand.get("vertebrate_like"), hmm_ok, at_locus=syn_ok)
    cand["name"] = base + cand.pop("name_suffix", "")
    cand["hmm_supported"] = hmm_ok
    cand["evidence_source"] = "+".join(label for label, ok in (("Synteny", syn_ok), ("HMM", hmm_ok)) if ok)
    # E1: in silico, supported by homology (MT-RNR2 position or Pfam GA)
    cand["evidence_tier"] = "E1" if syn_ok or hmm_ok else None


def run_pipeline(
    fasta_path: Path,
    output_prefix: str,
    table: int | str = "auto",
    enable_hmm: bool = True,
    enable_rescue: bool = True,
    enable_shlp: bool = True,
    circular: bool = False,
    all_candidates: bool = False,
    enable_synteny: bool = True,
) -> list[dict[str, Any]]:
    """End-to-end execution of the HumaninFinder pipeline.

    A Humanin frame is reported when it lies at a Humanin locus located by its MT-RNR2 position
    (synteny profiles) or carries a PF15040 domain at the gathering threshold (``_annotate_evidence``).
    """
    lint_report = lint_fasta(fasta_path, PIPELINE_LINT_OPTIONS)
    if not lint_report["summary"]["valid"]:
        errors = [i for i in lint_report["issues"] if i["severity"] == "error"]
        msg = f"FASTA lint failed for {fasta_path}: {len(errors)} errors."
        for err in errors[:3]:
            msg += f"\n  - [{err['code']}] {err['message']}"
        raise click.ClickException(msg)

    click.echo(f"[*] HumaninFinder v{__version__}: {' '.join(sys.argv)}")
    click.echo(f"[*] Processing sequences from {fasta_path}...")
    gathering = get_gathering_threshold()

    click.echo(f"\n{'-' * 101}")
    hdr = (
        f"{'Sequence ID':<22} | {'Peptide':<22} | {'Status':<12} | "
        f"{'1-Based Coords':<18} | {'Locus bits':<10} | {'HMM Score'}"
    )
    click.echo(hdr)
    click.echo(f"{'-' * 101}")

    requested_tables = parse_genetic_codes(table)
    multi_table = len(requested_tables) > 1

    all_accepted_results: list[dict[str, Any]] = []

    gap_rows: list[dict[str, Any]] = []

    with open_fasta_text(fasta_path) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            rec_id = str(record.id)
            full_seq = str(record.seq).upper()
            seq_len = len(full_seq)
            is_large_genome = seq_len > LARGE_SEQUENCE_BP
            # Only mitogenome-sized records are read as circular molecules
            record_circular = circular and not is_large_genome
            table_ids = requested_tables or [auto_genetic_code(record.description, seq_len)]

            if is_large_genome and not (enable_hmm or enable_synteny):
                click.echo(
                    f"{rec_id:<22} | [!] Skipped: sequences > {LARGE_SEQUENCE_BP} bp are scanned through the "
                    "synteny profile or the HMM-gated locus zoom; enable --synteny or --hmm."
                )
                continue

            candidates: list[dict[str, Any]] = []
            all_hmm_hits: list[dict[str, Any]] = []

            for t_id in table_ids:
                suffix = f"_T{t_id}" if multi_table else ""

                # 1. In-memory HMM search via chunked pyhmmer
                hmm_hits = []
                if enable_hmm:
                    hmm_hits = run_hmm_search_pyhmmer(full_seq, table=t_id)
                    for h in hmm_hits:
                        h["table"] = t_id
                    all_hmm_hits.extend(hmm_hits)

                # Humanin loci by their position in MT-RNR2 (multi-clade locus profile): the mitochondrial
                # locus and, in large records, every NUMT copy above the evidence threshold
                if enable_synteny:
                    candidates.extend(
                        _make_candidate(rec_id, locus, t_id, suffix, "Synteny")
                        for locus in locate_humanin_loci(full_seq, table=t_id, circular=record_circular)
                    )

                if is_large_genome:
                    # Two-stage HMM-gated locus zoom for nuclear chromosomes / large contigs. Only
                    # domains at the Pfam gathering threshold open a locus: sub-GA hits occur in random
                    # sequence (~1 per Mb).
                    gate_hits = [h for h in hmm_hits if h["score"] >= gathering]
                    if not gate_hits:
                        continue
                    candidate_regions = _merge_regions(
                        [
                            (max(0, h["dna_start"] - LOCUS_ZOOM_PAD_BP), min(seq_len, h["dna_end"] + LOCUS_ZOOM_PAD_BP))
                            for h in gate_hits
                        ]
                    )

                    raw_sorfs = []
                    for reg_start, reg_end in candidate_regions:
                        for s in find_sorfs(full_seq[reg_start:reg_end], table=t_id, circular=False):
                            s["start"] += reg_start
                            s["end"] += reg_start
                            raw_sorfs.append(s)
                    sorf_mode = "Scan_Zoom"
                else:
                    # Standard comprehensive scan for mitochondrial genomes (<= 50 kb)
                    raw_sorfs = find_sorfs(full_seq, table=t_id, circular=circular)
                    sorf_mode = "Scan"

                rescued_peps = []
                if enable_rescue and hmm_hits:
                    rescue_hits = hmm_hits
                    if is_large_genome:
                        rescue_hits = [
                            h
                            for h in hmm_hits
                            if any(lo <= h["dna_start"] and h["dna_end"] <= hi for lo, hi in candidate_regions)
                        ]
                    rescued_peps = guided_evolutionary_rescue(full_seq, rescue_hits, table=t_id)

                candidates.extend(_make_candidate(rec_id, s, t_id, suffix, sorf_mode) for s in raw_sorfs)
                candidates.extend(_make_candidate(rec_id, r, t_id, suffix, "Rescue") for r in rescued_peps)

            # 2. Link candidates with the best overlapping HMM domain
            if candidates:
                for cand in candidates:
                    c_start, c_end = cand["start_0based"], cand["end_0based"]
                    overlapping = [
                        h
                        for h in all_hmm_hits
                        if h["strand"] == cand["strand"]
                        and h.get("table") == cand["genetic_code"]
                        and calculate_overlap(c_start, c_end, h["dna_start"], h["dna_end"]) > 0
                    ]
                    if overlapping:
                        best_hmm = max(overlapping, key=lambda h: h["score"])
                        if cand.get("hmm_bitscore") is None or best_hmm["score"] > cand["hmm_bitscore"]:
                            cand["hmm_bitscore"] = float(best_hmm["score"])
                            cand["hmm_evalue"] = float(best_hmm["e_value"])
                    _annotate_evidence(cand, gathering)

            # 3. Non-redundant filtering
            nr_cands = filter_non_redundant(candidates)

            # 4. SHLP1–6 by their MT-RNR2 position, each in its primary (literature) translation context
            shlp_rows: list[dict[str, Any]] = []
            if enable_shlp and enable_synteny:
                # A start codon inside a deletion or an assembly gap yields no frame, so the peptide
                # is absent from the table. Record that, or it cannot be told from a missing locus.
                for gap in locus_gap_report(full_seq, circular=record_circular):
                    gap_rows.append({"sample_id": rec_id, **gap})
                # In a nuclear record SHLPs of mitochondrial primary context are read under the
                # nuclear code: a NUMT copy is not translated by the mitochondrial machinery. A record
                # is nuclear when it is large or when no requested code is mitochondrial (e.g. -g 1).
                nuclear_codes = [t for t in table_ids if t not in MITOCHONDRIAL_TABLES]
                is_nuclear = is_large_genome or len(nuclear_codes) == len(table_ids)
                nuclear_table = (nuclear_codes or [1])[0] if is_nuclear else None
                for sh in annotate_shlps(full_seq, circular=record_circular, nuclear_table=nuclear_table):
                    shlp_rows.append(
                        {
                            "sample_id": rec_id,
                            "name": sh["name"],
                            "status": sh["status"],
                            "seq": sh["seq"],
                            "length": sh["length"],
                            "strand": sh["strand"],
                            "frame": sh["frame"],
                            "start_0based": sh["start"],
                            "end_0based": sh["end"],
                            "start_1based": sh["start"] + 1,
                            "end_1based": sh["end"],
                            "circular_wrap": sh["circular_wrap"],
                            "hmm_bitscore": None,
                            "hmm_evalue": None,
                            "hmm_supported": False,
                            "synteny_score": sh["synteny_score"],
                            "locus_profile": sh["locus_profile"],
                            "locus_density": sh["locus_density"],
                            "start_offset_codons": sh["start_offset_codons"],
                            "evidence_source": "Synteny",
                            "evidence_tier": sh["evidence_tier"],
                            "reference_evidence_tier": sh["reference_evidence_tier"],
                            "identity_to_reference": sh["identity_to_reference"],
                            "genetic_code": sh["translation_table"],
                            "locus_mode": "Synteny",
                            "table_1_pep": "",
                        }
                    )

            reportable = [c for c in nr_cands if c["evidence_source"]]

            if not reportable and not shlp_rows:
                click.echo(f"{rec_id:<22} | >> No candidate sORFs detected.")
                continue

            # 5. Selection of candidates to report
            if all_candidates:
                selected = list(reportable)
            else:
                # Primary hit: best Humanin candidate per table
                selected = []
                for t_id in table_ids:
                    hn_cands_t = [c for c in reportable if c["genetic_code"] == t_id]
                    if hn_cands_t:
                        selected.append(
                            max(
                                hn_cands_t,
                                key=lambda x: (
                                    1.0 if x.get("synteny_score") is not None else 0.0,
                                    1.0 if x["hmm_supported"] else 0.0,
                                    STATUS_RANK.get(x["status"], 0.0),
                                    x.get("hmm_bitscore") or 0.0,
                                ),
                            )
                        )
            selected.extend(shlp_rows)

            for cand in selected:
                cand["seq_length"] = seq_len
                cand["seq_circular"] = record_circular
                coords = f"{cand['start_1based']}-{cand['end_1based']}"
                click.echo(
                    f"{rec_id:<22} | {cand.get('name', 'Humanin'):<22} | {cand['status']:<12} | "
                    f"{coords:<18} | {_fmt(cand.get('synteny_score'), '<10.0f', '-' * 10)} | "
                    f"{_fmt(cand.get('hmm_bitscore'), '.1f', '-')}"
                )
                all_accepted_results.append(cand)

    write_deterministic_outputs(all_accepted_results, output_prefix)
    if gap_rows:
        gap_cols = ["sample_id", "name", "locus_profile", "synteny_score", "start_column", "gap_columns"]
        gap_sorted = sorted(gap_rows, key=lambda r: (r["sample_id"], r["name"], r["locus_profile"]))
        write_csv(gap_sorted, gap_cols, f"{output_prefix}_gaps.csv")
    write_run_provenance(
        output_prefix,
        fasta_path,
        {
            "table": table,
            "hmm": enable_hmm,
            "rescue": enable_rescue,
            "shlp": enable_shlp,
            "synteny": enable_synteny,
            "circular": circular,
            "all_candidates": all_candidates,
        },
        len(all_accepted_results),
        gap_rows=len(gap_rows),
    )
    click.echo(f"\n[+] HumaninFinder v{__version__} Execution Complete:")
    click.echo(f"    - CSV Table:   {output_prefix}_results.csv")
    click.echo(f"    - FASTA Track: {output_prefix}_results.fasta")
    click.echo(f"    - GFF3 Track:  {output_prefix}_results.gff3")
    click.echo(f"    - BED Track:   {output_prefix}_results.bed")
    click.echo(f"    - Run record:  {output_prefix}_run.json")
    if gap_rows:
        affected = sorted({r["sample_id"] for r in gap_rows})
        peptides = sorted({r["name"] for r in gap_rows})
        click.echo(f"    - Start-codon gaps: {output_prefix}_gaps.csv")
        click.echo(
            f"[!] {', '.join(peptides)} absent from {len(affected)} record(s): the start codon falls in a "
            "gap, so no frame could be read. Absence here is not evidence that the locus is missing."
        )

    return all_accepted_results


def _pipeline_options(func):
    """Options shared by the predict and scan commands."""
    options = [
        click.option(
            "--input",
            "-i",
            required=True,
            type=click.Path(exists=True, dir_okay=False, path_type=Path),
            help="Input FASTA file (plain or gzip-compressed).",
        ),
        click.option("--output", "-o", required=True, help="Output prefix path."),
        click.option(
            "--table",
            "-g",
            type=str,
            default="auto",
            help="NCBI Genetic Code Table(s), e.g. 2 or 2,1 (comma-separated), or 'auto': Table 2 for "
            "mitochondrial records (header or <= 50 kb), Table 1 otherwise (vertebrate codes). Default: auto.",
        ),
        click.option("--hmm/--no-hmm", default=True, help="Enable/disable pyhmmer HMM profile search."),
        click.option("--rescue/--no-rescue", default=True, help="Enable/disable guided evolutionary rescue."),
        click.option("--shlp/--no-shlp", default=True, help="Enable/disable SHLP1-6 micropeptide annotation."),
        click.option(
            "--synteny/--no-synteny",
            default=True,
            help="Enable/disable location of Humanin loci by their MT-RNR2 position (multi-clade locus profile).",
        ),
        click.option("--circular", is_flag=True, help="Enable circular chromosome wrap-around search."),
    ]
    for option in reversed(options):
        func = option(func)
    return func


@click.group()
@click.version_option(version=__version__)
def main():
    """HumaninFinder: homology-based location and annotation of Humanin and SHLP1-6."""
    pass


@main.command()
@_pipeline_options
@click.option(
    "--all-candidates", is_flag=True, help="Report every non-redundant candidate instead of the best one per table."
)
def predict(
    input,
    output,
    table,
    hmm,
    rescue,
    shlp,
    synteny,
    circular,
    all_candidates,
):
    """Predict and annotate Humanin and SHLPs in target mitochondrial sequences (best hit per table)."""
    run_pipeline(
        fasta_path=input,
        output_prefix=output,
        table=table,
        enable_hmm=hmm,
        enable_rescue=rescue,
        enable_shlp=shlp,
        circular=circular,
        all_candidates=all_candidates,
        enable_synteny=synteny,
    )


@main.command()
@_pipeline_options
@click.option(
    "--all-candidates/--best-only",
    default=True,
    help="Report every non-redundant candidate (default) or only the best one per table.",
)
def scan(
    input,
    output,
    table,
    hmm,
    rescue,
    shlp,
    synteny,
    circular,
    all_candidates,
):
    """Scan whole mitogenomes, contigs, or nuclear genomes (NUMTs), reporting every locus."""
    run_pipeline(
        fasta_path=input,
        output_prefix=output,
        table=table,
        enable_hmm=hmm,
        enable_rescue=rescue,
        enable_shlp=shlp,
        circular=circular,
        all_candidates=all_candidates,
        enable_synteny=synteny,
    )


@main.command()
@click.option("--input", "-i", required=True, type=click.Path(exists=True, path_type=Path), help="Cohort FASTA.")
@click.option("--output", "-o", required=True, help="Output prefix.")
@click.option(
    "--population",
    is_flag=True,
    help="Cohort of one species: count segregating sites and align globally by rCRS coordinates.",
)
@click.option("--workers", default=1, show_default=True, type=click.IntRange(1, 64), help="Parallel processes.")
@click.option("--limit", default=0, type=click.IntRange(0), help="Use only the first N sequences (0 = all).")
@click.option(
    "--reference",
    type=click.Path(exists=True, path_type=Path),
    default=None,
    help="Override the bundled rCRS reference genes.",
)
def conservation(input, output, population, workers, limit, reference):
    """Test whether conservation at the Humanin window tracks the peptide or the host 16S rRNA.

    Runs the codon-position, shifted-frame and rRNA-null controls over a cohort of genomes, with
    MT-CO1 as the positive control. This command does NOT use the MT-RNR2 profiles or the Pfam
    domain: it aligns each genome to the bundled rCRS reference genes, so the result is independent
    of how the locus caller behaves.
    """
    click.echo(f"[*] HumaninFinder v{__version__}: {' '.join(sys.argv)}")
    click.echo(f"[*] Coding-evidence controls over {input}...")
    record = analyse_cohort(
        input,
        population=population,
        reference_path=reference,
        workers=workers,
        limit=limit,
    )
    path = write_report(record, output)
    click.echo("")
    click.echo(format_summary(record))
    write_run_provenance(
        output,
        input,
        {"population": population, "workers": workers, "limit": limit, "reference": str(reference or "bundled")},
        record["genomes"],
    )
    click.echo("")
    click.echo(f"[+] Report:     {path}")
    click.echo(f"[+] Run record: {output}_run.json")


@main.command()
@click.argument("input", type=click.Path(exists=True, path_type=Path))
def lint(input):
    """Validate FASTA formatting and compliance without modifying it."""
    report = lint_fasta(input, PIPELINE_LINT_OPTIONS)
    click.echo(f"[*] Linting {input}...")
    is_valid = report["summary"]["valid"]
    record_count = report["summary"]["record_count"]
    if is_valid:
        click.secho(f"[+] FASTA is 100% valid ({record_count} records).", fg="green")
    else:
        errors = [i for i in report["issues"] if i["severity"] == "error"]
        warnings = [i for i in report["issues"] if i["severity"] == "warning"]
        click.secho(f"[-] FASTA has {len(errors)} errors and {len(warnings)} warnings:", fg="red")
        for err in errors:
            click.echo(f"    - Error: [{err['code']}] {err['message']}")
        for warn in warnings:
            click.echo(f"    - Warning: [{warn['code']}] {warn['message']}")
        sys.exit(1)


@main.command()
def setup():
    """Verify system prerequisites and every bundled data file.

    Checks all of them, not only the locus profiles: a package whose PF15040.hmm, SHLP references or
    conservation reference failed to install would otherwise pass this command and then fail, or
    silently under-report, in the middle of a run.
    """
    click.echo(f"[*] Verifying HumaninFinder v{__version__} prerequisites...")
    import pyhmmer

    from .conservation import DEFAULT_REFERENCE, GENE_COORDS, load_reference
    from .core import DATA_DIR, DEFAULT_HUMANIN_HMM, get_gathering_threshold, load_locus_profiles
    from .shlp import SHLP_REF_PATH, load_shlp_references

    click.secho(f"[+] pyhmmer (in-memory HMMER {pyhmmer.__version__}) verified.", fg="green")

    failures: list[str] = []

    # Every file the package ships, and what each one must contain to be usable
    if not DEFAULT_HUMANIN_HMM.exists():
        failures.append(f"missing {DEFAULT_HUMANIN_HMM.name}")
    else:
        ga = get_gathering_threshold()
        click.secho(f"[+] PF15040 Humanin domain: gathering threshold {ga:.1f} bits.", fg="green")

    if not SHLP_REF_PATH.exists():
        failures.append(f"missing {SHLP_REF_PATH.name}")
    else:
        shlps = load_shlp_references().get("shlp_peptides", {})
        if len(shlps) != 6:
            failures.append(f"{SHLP_REF_PATH.name} holds {len(shlps)} peptides, expected 6")
        else:
            click.secho(f"[+] SHLP references: {', '.join(sorted(shlps))}.", fg="green")

    try:
        ref = load_reference()
        missing = sorted(set(GENE_COORDS) - set(ref))
        if missing:
            failures.append(f"conservation reference lacks {', '.join(missing)}")
        else:
            sizes = ", ".join(f"{g} {len(ref[g])} nt" for g in sorted(ref))
            click.secho(f"[+] Conservation reference ({DEFAULT_REFERENCE.name}): {sizes}.", fg="green")
    except (OSError, ValueError) as err:
        failures.append(f"conservation reference unreadable: {err}")

    orphans = sorted(
        p.name
        for p in DATA_DIR.iterdir()
        if p.is_file()
        and p.name
        not in {
            DEFAULT_HUMANIN_HMM.name,
            SHLP_REF_PATH.name,
            DEFAULT_REFERENCE.name,
            "MT-RNR2.hmm",
            "Humanin_core.hmm",
            "mdp_locus_profiles.json",
        }
    )
    if orphans:
        click.secho(f"[!] Unexpected files in data/: {', '.join(orphans)}", fg="yellow")

    profiles = load_locus_profiles()
    click.secho(
        "[+] MT-RNR2 locus profiles: "
        + ", ".join(
            f"{n} ({p['profile_length']} states, >= {p['evidence_bitscore']:.0f} bits)" for n, p in profiles.items()
        ),
        fg="green",
    )
    click.secho(
        "[+] Vertebrate-like density (Humanin vs Humanin_position): "
        + ", ".join(
            f"{n} >= {p['vertebrate_min_density']:.3f} bits/column"
            for n, p in profiles.items()
            if p.get("vertebrate_min_density") is not None
        ),
        fg="green",
    )
    if failures:
        for f in failures:
            click.secho(f"[-] {f}", fg="red")
        raise click.ClickException(
            f"{len(failures)} bundled data file(s) missing or unusable. Reinstall the package; "
            "a partial install gives wrong results rather than an error."
        )
    click.secho("\n[!] Setup inspection complete: all bundled data files verified.", fg="green", bold=True)


if __name__ == "__main__":
    main()
