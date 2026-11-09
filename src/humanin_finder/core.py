"""Core genomic algorithms, sORF discovery, 16S locus detection, and pyhmmer in-memory search."""

from __future__ import annotations

import json
import logging
from functools import lru_cache
from pathlib import Path
from typing import Any

import pyhmmer
from Bio.Data import CodonTable
from Bio.Seq import Seq

logger = logging.getLogger(__name__)

DATA_DIR = Path(__file__).resolve().parent / "data"
DEFAULT_HUMANIN_HMM = DATA_DIR / "PF15040.hmm"
# Multi-clade vertebrate DNA profiles of MT-RNR2 (whole 16S and Humanin core; built by the
# reproducibility repository, 10_build_locus_profile.py). The .json holds, per profile, the match
# state of each peptide start codon, the human reference ORFs and the calibrated bit thresholds.
DEFAULT_LOCUS_PROFILES = DATA_DIR / "mdp_locus_profiles.json"
# A stop codon before this fraction of the human peptide length interrupts the homologous frame
INTERRUPTION_FRACTION = 0.75
# The species' own ORF may start up to this many codons from the mapped human start codon (both
# directions), so its status is judged on the local reading frame, not on human coordinates
LOCAL_START_WINDOW = 6
# NCBI genetic codes of mitochondria (all start codons of the table are accepted as initiators)
MITOCHONDRIAL_TABLES = frozenset({2, 3, 4, 5, 9, 13, 14, 16, 21, 22, 23, 24, 33})
# Pfam gathering threshold (GA) of PF15040, used when the HMM file does not declare one
FALLBACK_GATHERING_THRESHOLD = 27.0

# Structural status of a reading frame (no homology or functional claim is implied)
STATUS_COMPLETE = "Complete"  # initiator Met through an in-frame stop codon
STATUS_INCOMPLETE = "Incomplete"  # no initiator Met and/or no terminating stop inside the window
STATUS_INTERRUPTED = "Interrupted"  # internal in-frame stop codon under the evaluated genetic code


@lru_cache(maxsize=8)
def _load_hmm(hmm_path: str) -> pyhmmer.plan7.HMM:
    """Read (once per process) the first profile HMM stored in a file."""
    with pyhmmer.plan7.HMMFile(hmm_path) as hmm_file:
        return hmm_file.read()


def get_gathering_threshold(hmm_path: str | Path | None = None) -> float:
    """Return the per-domain Pfam gathering threshold (GA, bits) of the HMM profile."""
    hmm_path = Path(hmm_path) if hmm_path is not None else DEFAULT_HUMANIN_HMM
    if not hmm_path.exists():
        return FALLBACK_GATHERING_THRESHOLD
    hmm = _load_hmm(str(hmm_path))
    if hmm.cutoffs.gathering_available():
        return float(hmm.cutoffs.gathering[1])
    return FALLBACK_GATHERING_THRESHOLD


def get_translation_table(table: int | str = 2) -> CodonTable.NCBICodonTableDNA:
    """Retrieve an NCBI DNA codon table by integer ID (1-33)."""
    try:
        table_id = int(table)
        return CodonTable.unambiguous_dna_by_id[table_id]
    except (KeyError, ValueError) as err:
        raise ValueError(
            f"Invalid NCBI translation table ID: {table}. Not every integer in 1-33 is defined "
            "(7, 8, 15 and 17-20 are vacant); mitochondrial tables are "
            f"{sorted(MITOCHONDRIAL_TABLES)}."
        ) from err


def translate_dna(
    dna_seq: str | Seq,
    table: int | str = 2,
    to_stop: bool = False,
) -> str:
    """Translate DNA to amino acid sequence with safe trimming of partial terminal codons."""
    seq_str = str(dna_seq).upper()
    if not seq_str:
        return ""

    trim_len = (len(seq_str) // 3) * 3
    if trim_len == 0:
        return ""

    seq_obj = Seq(seq_str[:trim_len])
    table_id = int(table)
    return str(seq_obj.translate(table=table_id, to_stop=to_stop))


def is_degenerate_peptide(seq: str, max_internal_stops: int = 2) -> bool:
    """True when a translated frame is not worth reporting as a candidate sORF.

    Four independent rejections: empty sequence; more than ``max_internal_stops`` stop codons
    (``*``); more than 25% ``X``; one residue above 60% of a sequence of at least 10 aa. Only the
    last is a compositional (low-complexity) test, which is why the former name ``is_low_complexity``
    was misleading.
    The stop rule also drops degraded NUMT frames, so it is applied only on the sORF route
    (``find_sorfs``), never to the MDP locus route.
    """
    clean_seq = seq.upper()
    if not clean_seq:
        return True

    if clean_seq.count("*") > max_internal_stops:
        return True

    if clean_seq.count("X") / len(clean_seq) > 0.25:
        return True

    counts: dict[str, int] = {}
    for aa in clean_seq:
        counts[aa] = counts.get(aa, 0) + 1
    max_freq = max(counts.values()) / len(clean_seq)
    if max_freq > 0.60 and len(clean_seq) >= 10:
        return True

    return False


def calculate_overlap(start1: int, end1: int, start2: int, end2: int) -> int:
    """Calculate overlap base pairs between two 0-based half-open intervals [start, end)."""
    overlap = min(end1, end2) - max(start1, start2)
    return max(0, overlap)


def find_sorfs(
    dna_seq: str | Seq,
    table: int | str = 2,
    min_len: int = 10,
    max_len: int = 60,
    allow_alternative_starts: bool = False,
    circular: bool = False,
    strand_filter: int | None = None,
) -> list[dict[str, Any]]:
    """Discover small open reading frames (sORFs) across reading frames.

    Coordinates are 0-based half-open intervals [start, end) on the forward strand.
    For ORFs crossing the origin of a circular sequence (``circular_wrap=True``),
    ``start`` lies in [0, len) and ``end = start + span`` exceeds the sequence length,
    following the GFF3 convention for circular landmarks.
    """
    seq_str = str(dna_seq).upper()
    seq_len = len(seq_str)
    if seq_len < min_len * 3:
        return []

    codon_tbl = get_translation_table(table)
    table_id = int(table)
    stop_codons = set(codon_tbl.stop_codons)
    standard_starts = set(get_translation_table(1).start_codons)

    if allow_alternative_starts:
        start_codons = set(codon_tbl.start_codons)
    else:
        start_codons = {c for c in codon_tbl.start_codons if translate_dna(c, table=table_id) == "M"}
        if not start_codons:
            start_codons = {"ATG"}

    sorfs = []
    strands: list[tuple[int, str]] = []
    if strand_filter is None or strand_filter == 1:
        strands.append((1, seq_str))
    if strand_filter is None or strand_filter == -1:
        rev_comp = str(Seq(seq_str).reverse_complement())
        strands.append((-1, rev_comp))

    for strand, nuc_strand in strands:
        search_seq = nuc_strand + (nuc_strand[: max_len * 3] if circular else "")
        strand_len = len(nuc_strand)

        for frame in range(3):
            frame_codons = [search_seq[i : i + 3] for i in range(frame, len(search_seq) - 2, 3)]
            num_codons = len(frame_codons)

            i = 0
            while i < num_codons:
                codon = frame_codons[i]
                start_nuc_pos = frame + i * 3

                if circular and start_nuc_pos >= strand_len:
                    break

                if codon in start_codons:
                    j = i + 1
                    found_stop = False
                    # ORFs longer than max_len are discarded, so the stop search is bounded
                    last_codon = min(num_codons, i + max_len + 1)
                    while j < last_codon:
                        if frame_codons[j] in stop_codons:
                            found_stop = True
                            break
                        j += 1

                    if found_stop:
                        pep_len = j - i
                        if min_len <= pep_len <= max_len:
                            end_nuc_pos = frame + (j + 1) * 3
                            coding_dna = search_seq[start_nuc_pos : frame + j * 3]
                            # Any start codon is decoded as the initiator methionine
                            pep_seq = "M" + translate_dna(coding_dna, table=table_id)[1:]

                            is_wrap = circular and (end_nuc_pos > strand_len)

                            if strand == 1:
                                d_start = start_nuc_pos
                                d_end = end_nuc_pos
                            elif is_wrap:
                                d_start = 2 * strand_len - end_nuc_pos
                                d_end = d_start + (end_nuc_pos - start_nuc_pos)
                            else:
                                d_start = strand_len - end_nuc_pos
                                d_end = strand_len - start_nuc_pos

                            # The cytosolic (Table 1) reading only exists if the codon also initiates there
                            dual_translations = {}
                            if table_id == 2 and codon in standard_starts:
                                dual_translations["table_1_pep"] = translate_dna(
                                    search_seq[start_nuc_pos : start_nuc_pos + (max_len + 1) * 3],
                                    table=1,
                                    to_stop=True,
                                )

                            sorfs.append(
                                {
                                    "seq": pep_seq,
                                    "dna_seq": coding_dna,
                                    "start": d_start,
                                    "end": d_end,
                                    "strand": strand,
                                    "frame": frame,
                                    "status": STATUS_COMPLETE,
                                    "length": len(pep_seq),
                                    "start_codon": codon,
                                    "stop_codon": frame_codons[j],
                                    "circular_wrap": is_wrap,
                                    "table": table_id,
                                    **dual_translations,
                                }
                            )
                i += 1

    return sorfs


def _search_chunk_pyhmmer(
    seq_str: str,
    hmm: pyhmmer.plan7.HMM,
    alphabet: pyhmmer.easel.Alphabet,
    table_id: int,
    evalue_threshold: float,
    bitscore_threshold: float,
) -> list[dict[str, Any]]:
    """Search a discrete nucleotide chunk across all 6 frames using pyhmmer."""
    seq_len = len(seq_str)
    if seq_len < 30:
        return []

    strands = [
        (1, seq_str),
        (-1, str(Seq(seq_str).reverse_complement())),
    ]

    easel_seqs: list[pyhmmer.easel.DigitalSequence] = []
    query_meta: list[dict[str, Any]] = []

    for strand, nuc in strands:
        for frame in range(3):
            frame_dna = nuc[frame:]
            pep = translate_dna(frame_dna, table=table_id)
            if len(pep) < 10:
                continue

            # 1. Full-frame query
            q_name = f"s{strand}_f{frame}_full"
            easel_seqs.append(
                pyhmmer.easel.TextSequence(name=q_name.encode("utf-8"), sequence=pep).digitize(alphabet)
            )
            query_meta.append(
                {"strand": strand, "frame": frame, "offset": 0, "pep": pep, "name": q_name}
            )

            # 2. Overlapping 35-aa windows (stride 12) to detect pseudogenes with internal stops
            if len(pep) > 30:
                for w_idx in range(0, len(pep) - 15, 12):
                    chunk = pep[w_idx : min(len(pep), w_idx + 35)]
                    q_name = f"s{strand}_f{frame}_w{w_idx}"
                    easel_seqs.append(
                        pyhmmer.easel.TextSequence(name=q_name.encode("utf-8"), sequence=chunk).digitize(alphabet)
                    )
                    query_meta.append(
                        {"strand": strand, "frame": frame, "offset": w_idx, "pep": pep, "name": q_name}
                    )

    if not easel_seqs:
        return []

    search_results = list(pyhmmer.hmmer.hmmsearch(hmm, easel_seqs, E=evalue_threshold, cpus=1))
    if not search_results:
        return []

    meta_by_name = {m["name"]: m for m in query_meta}
    raw_hits: list[dict[str, Any]] = []

    for batch in search_results:
        for hit in batch:
            if hit.score < bitscore_threshold:
                continue
            h_name = hit.name.decode("utf-8") if isinstance(hit.name, bytes) else str(hit.name)
            meta_info = meta_by_name.get(h_name)
            if not meta_info:
                continue

            strand = meta_info["strand"]
            frame = meta_info["frame"]
            offset = meta_info["offset"]
            pep = meta_info["pep"]

            for domain in hit.domains:
                if domain.score < bitscore_threshold:
                    continue

                ali_from = offset + domain.alignment.target_from - 1  # 0-based target peptide start
                ali_to = offset + domain.alignment.target_to  # 0-based target peptide end

                nuc_start = ali_from * 3 + frame
                nuc_end = ali_to * 3 + frame

                if strand == 1:
                    d_start = nuc_start
                    d_end = nuc_end
                else:
                    d_start = seq_len - nuc_end
                    d_end = seq_len - nuc_start

                matched_pep = pep[ali_from:ali_to]

                raw_hits.append(
                    {
                        "strand": strand,
                        "frame": frame,
                        "score": float(domain.score),
                        "e_value": float(domain.i_evalue),
                        "c_evalue": float(domain.c_evalue),
                        "aa_start": ali_from,
                        "aa_end": ali_to,
                        "dna_start": d_start,
                        "dna_end": d_end,
                        "matched_pep": matched_pep,
                    }
                )

    return raw_hits


def run_hmm_search_pyhmmer(
    dna_seq: str | Seq,
    hmm_path: str | Path | None = None,
    table: int | str = 2,
    evalue_threshold: float = 1.0,
    bitscore_threshold: float = 5.0,
    chunk_size: int = 60_000,
    chunk_overlap: int = 3_000,
) -> list[dict[str, Any]]:
    """Perform in-memory HMMER search across all 6 reading frames with constant O(1) memory chunking."""
    if hmm_path is None:
        hmm_path = DEFAULT_HUMANIN_HMM

    hmm_path = Path(hmm_path)
    if not hmm_path.exists():
        logger.warning("HMM file not found at: %s", hmm_path)
        return []

    seq_str = str(dna_seq).upper()
    seq_len = len(seq_str)
    if seq_len < 30:
        return []

    table_id = int(table)
    alphabet = pyhmmer.easel.Alphabet.amino()
    hmm = _load_hmm(str(hmm_path))

    raw_hits: list[dict[str, Any]] = []

    if seq_len <= chunk_size:
        raw_hits = _search_chunk_pyhmmer(
            seq_str, hmm, alphabet, table_id, evalue_threshold, bitscore_threshold
        )
    else:
        step = chunk_size - chunk_overlap
        for c_start in range(0, seq_len, step):
            c_end = min(seq_len, c_start + chunk_size)
            chunk_nuc = seq_str[c_start:c_end]
            if len(chunk_nuc) < 30:
                continue
            c_hits = _search_chunk_pyhmmer(
                chunk_nuc, hmm, alphabet, table_id, evalue_threshold, bitscore_threshold
            )
            for h in c_hits:
                h["dna_start"] += c_start
                h["dna_end"] += c_start
                raw_hits.append(h)
            if c_end == seq_len:
                break

    # Deduplicate overlapping hits on same strand keeping highest score
    raw_hits.sort(key=lambda h: h["score"], reverse=True)
    hits: list[dict[str, Any]] = []
    for h in raw_hits:
        is_dup = False
        h_len = max(1, h["dna_end"] - h["dna_start"])
        for kept in hits:
            if kept["strand"] == h["strand"]:
                overlap = calculate_overlap(h["dna_start"], h["dna_end"], kept["dna_start"], kept["dna_end"])
                if overlap / h_len > 0.50:
                    is_dup = True
                    break
        if not is_dup:
            hits.append(h)

    return hits


def guided_evolutionary_rescue(
    dna_seq: str | Seq,
    hmm_hits: list[dict[str, Any]],
    table: int | str = 2,
    min_len: int = 15,
    max_len: int = 45,
    max_internal_stops: int = 2,
) -> list[dict[str, Any]]:
    """Rescue incomplete or interrupted Humanin-like reading frames guided by HMM evidence.

    An internal stop codon marks the frame as ``Interrupted`` under the evaluated genetic code; it
    is not by itself evidence of pseudogenisation (the peptide may be produced under another code).
    """
    seq_str = str(dna_seq).upper()
    seq_len = len(seq_str)
    if not hmm_hits or seq_len < min_len * 3:
        return []

    table_id = int(table)
    rescued = []
    seen_peptides = set()

    for hit in hmm_hits:
        strand = hit["strand"]
        pad = 300
        w_start = max(0, hit["dna_start"] - pad)
        w_end = min(seq_len, hit["dna_end"] + pad)
        window_dna = seq_str[w_start:w_end]

        if strand == 1:
            nuc_source = window_dna
            target_start = hit["dna_start"] - w_start
            target_end = hit["dna_end"] - w_start
        else:
            nuc_source = str(Seq(window_dna).reverse_complement())
            target_start = w_end - hit["dna_end"]
            target_end = w_end - hit["dna_start"]

        frame = target_start % 3
        frame_dna = nuc_source[frame:]
        full_aa = translate_dna(frame_dna, table=table_id)
        aa_start = (target_start - frame) // 3
        aa_end = (target_end - frame) // 3

        best_start_aa = aa_start
        for k in range(aa_start, max(-1, aa_start - 15), -1):
            if k >= 0 and full_aa[k] == "M":
                if full_aa[k:aa_start].count("*") <= max_internal_stops:
                    best_start_aa = k
                    break

        first_stop = full_aa.find("*", best_start_aa)
        if first_stop != -1 and (first_stop - best_start_aa) >= min_len:
            best_end_aa = min(best_start_aa + max_len, first_stop)
        elif first_stop != -1 and (first_stop - best_start_aa) < min_len:
            # Internal stop under this code (interrupted frame); capture the ~21-25 aa domain
            second_stop = full_aa.find("*", first_stop + 1)
            target_span = best_start_aa + 21
            if second_stop != -1 and min_len <= (second_stop - best_start_aa) <= max_len:
                best_end_aa = second_stop
            else:
                best_end_aa = min(len(full_aa), max(target_span, aa_end))
        else:
            best_end_aa = min(len(full_aa), best_start_aa + 24)

        cand_pep = full_aa[best_start_aa:best_end_aa]

        if len(cand_pep) < min_len or len(cand_pep) > max_len:
            center = (aa_start + aa_end) // 2
            half = min_len // 2
            w_s = max(0, center - half)
            w_e = min(len(full_aa), w_s + min_len)
            cand_pep = full_aa[w_s:w_e]
            best_start_aa = w_s
            best_end_aa = w_e

        if len(cand_pep) < min_len or is_degenerate_peptide(cand_pep, max_internal_stops):
            continue

        if cand_pep in seen_peptides:
            continue
        seen_peptides.add(cand_pep)

        if "*" in cand_pep:
            status = STATUS_INTERRUPTED
        elif cand_pep.startswith("M") and first_stop == best_end_aa:
            status = STATUS_COMPLETE
        else:
            status = STATUS_INCOMPLETE

        nuc_start = best_start_aa * 3 + frame
        nuc_end = best_end_aa * 3 + frame

        if strand == 1:
            d_start = w_start + nuc_start
            d_end = w_start + nuc_end
        else:
            d_start = w_end - nuc_end
            d_end = w_end - nuc_start

        overlap_bp = calculate_overlap(d_start, d_end, hit["dna_start"], hit["dna_end"])

        rescued.append(
            {
                "seq": cand_pep,
                "start": d_start,
                "end": d_end,
                "strand": strand,
                "frame": frame,
                "status": status,
                "length": len(cand_pep),
                "hmm_score": hit["score"],
                "hmm_evalue": hit["e_value"],
                "hmm_overlap_bp": overlap_bp,
                "table": table_id,
                "rescued": True,
            }
        )

    return rescued


def score_peptide_hmm(peptide: str, hmm_path: str | Path | None = None) -> float:
    """Best PF15040 domain bitscore of a single peptide (stop codons read as unknown residues)."""
    hmm_path = Path(hmm_path) if hmm_path is not None else DEFAULT_HUMANIN_HMM
    pep = peptide.replace("*", "X")
    if len(pep) < 5 or not hmm_path.exists():
        return 0.0
    alphabet = pyhmmer.easel.Alphabet.amino()
    seq = pyhmmer.easel.TextSequence(name=b"locus", sequence=pep).digitize(alphabet)
    hits = pyhmmer.hmmer.hmmsearch(_load_hmm(str(hmm_path)), [seq], E=1000.0, Z=1, cpus=1)
    return max((domain.score for top in hits for hit in top for domain in hit.domains), default=0.0)


def _read_locus_frame(
    oriented: str,
    t_start: int,
    strand: int,
    seq_len: int,
    score: float,
    reference_orf: str,
    table: int | str,
    max_len: int,
    score_pf15040: bool = True,
) -> dict[str, Any] | None:
    """Read the species' own reading frame at a mapped peptide start codon (see ``locate_mdp_loci``).

    ``reference_orf`` is the human ORF (start codon to last sense codon) that defines the homologous
    span; it only sets the span length and the descriptive ``identity_to_reference``.
    """
    table_id = int(table)
    codon_tbl = get_translation_table(table_id)
    stops = set(codon_tbl.stop_codons)
    starts = set(codon_tbl.start_codons)
    if table_id not in MITOCHONDRIAL_TABLES:
        # Outside organelles only Met-decoding initiators are accepted (no near-cognate CTG/TTG starts)
        starts = {c for c in starts if translate_dna(c, table=table_id) == "M"}
    human_len = len(reference_orf) // 3
    codons = [oriented[i : i + 3] for i in range(t_start, min(len(oriented) - 2, t_start + 3 * (max_len + 1)), 3)]
    if len(codons) < human_len // 2:
        return None
    # In-frame codons upstream of the mapped start, nearest first, up to the first stop codon
    upstream = []
    for k in range(1, LOCAL_START_WINDOW + 1):
        pos = t_start - 3 * k
        if pos < 0 or oriented[pos : pos + 3] in stops:
            break
        upstream.append(oriented[pos : pos + 3])

    def local_orf(k: int) -> tuple[int, int | None] | None:
        """ORF starting k codons from the mapped start: (k, stop index) if it covers enough of the core."""
        first = max(k, 0)
        # An upstream start also reads through the mapped codon, which may itself be a stop
        stop = next((j for j in range(first + (k >= 0), len(codons)) if codons[j] in stops), None)
        if stop is not None and stop - k > max_len:
            return None
        covered = min(stop if stop is not None else len(codons), human_len) - first
        return (k, stop) if covered >= INTERRUPTION_FRACTION * human_len else None

    # Nearest start codon to the mapped one; an upstream start wins ties (it keeps the whole core)
    offsets = [0] + [k for d in range(1, LOCAL_START_WINDOW + 1) for k in (-d, d)]
    for k in offsets:
        if k < 0:
            codon = upstream[-k - 1] if -k <= len(upstream) else None
        else:
            codon = codons[k] if k < len(codons) else None
        if codon in starts and (orf := local_orf(k)) is not None:
            break
    else:
        orf = None
    # Homologous span from the mapped start (the human-frame reading, reported when no local ORF exists)
    stop_idx = next((i for i, c in enumerate(codons) if c in stops), None)

    if orf is not None:
        k, stop = orf
        orf_codons = (upstream[:-k][::-1] if k < 0 else []) + codons[
            max(k, 0) : stop if stop is not None else len(codons)
        ]
        peptide = "M" + translate_dna("".join(orf_codons), table=table_id)[1:]
        status = STATUS_COMPLETE if stop is not None else STATUS_INCOMPLETE
        span_start = 3 * k
        span_codons = len(orf_codons) + (1 if stop is not None else 0)
        start_codon = orf_codons[0]
    elif stop_idx is not None and stop_idx < INTERRUPTION_FRACTION * human_len:
        k, status = None, STATUS_INTERRUPTED
        span_codons = min(len(codons), human_len)
        peptide = translate_dna("".join(codons[:span_codons]), table=table_id)
        span_start, start_codon = 0, codons[0]
    else:
        k, status = None, STATUS_INCOMPLETE
        n_codons = stop_idx if stop_idx is not None else len(codons)
        peptide = translate_dna("".join(codons[:n_codons]), table=table_id)
        span_codons = n_codons + (1 if stop_idx is not None else 0)
        span_start, start_codon = 0, codons[0]

    # Comparison with the human reference is reported apart from the status: ungapped identity of
    # the homologous span (mapped start, human length) to the human peptide under the same code
    human_pep = translate_dna(reference_orf, table=table_id)
    target_pep = translate_dna("".join(codons[:human_len]), table=table_id)
    identity = sum(a == b for a, b in zip(human_pep, target_pep, strict=False)) / human_len

    o_start = int(t_start) + span_start
    o_end = o_start + 3 * span_codons
    if strand == 1:
        start, end = o_start, o_end
    else:
        start, end = seq_len - o_end, seq_len - o_start

    return {
        "seq": peptide,
        "start": start,
        "end": end,
        "strand": strand,
        "frame": o_start % 3,
        "status": status,
        "length": len(peptide),
        "start_codon": start_codon,
        "table": table_id,
        "synteny_score": float(score),
        "start_offset_codons": k,
        "identity_to_reference": round(identity, 4),
        "hmm_score": score_peptide_hmm(peptide) if score_pf15040 else None,
        "hmm_evalue": None,
        "rescued": True,
    }


@lru_cache(maxsize=2)
def load_locus_profiles(profiles_path: str | None = None) -> dict[str, dict[str, Any]]:
    """Read the MT-RNR2 locus profiles: {name: {"hmm": HMM, **metadata}} in search order."""
    path = Path(profiles_path) if profiles_path else DEFAULT_LOCUS_PROFILES
    meta = json.loads(path.read_text())
    profiles = {}
    for name, prof in meta["profiles"].items():
        with pyhmmer.plan7.HMMFile(str(path.parent / prof["file"])) as fh:
            profiles[name] = {**prof, "hmm": fh.read()}
    return profiles


def _profile_positions(
    seq_str: str, hmm: pyhmmer.plan7.HMM, states: dict[str, tuple[int, int]], threshold: float, cpus: int
) -> list[tuple[float, str, int, int, float]]:
    """nhmmer hits of a locus profile mapped to peptide start codons.

    ``states`` maps a peptide to (match state of the first base of its start codon, orientation of
    the peptide relative to the profile). Returns (bitscore, peptide, peptide strand on the target,
    0-based start on that strand, bits per aligned profile column) for every hit at or above
    ``threshold``.

    A start-codon column that aligns to a gap in the target yields no position, because the codon
    itself is not in the record. ``locus_gap_report`` reports those cases separately, so that a
    peptide missing from the output can be told apart from a locus the genome does not have.
    """
    alphabet = pyhmmer.easel.Alphabet.dna()
    target = pyhmmer.easel.TextSequence(name=b"target", sequence=seq_str).digitize(alphabet)
    seq_len = len(seq_str)
    by_state = {state: [] for state, _ in states.values()}
    for name, (state, orientation) in states.items():
        by_state[state].append((name, orientation))
    positions = []
    for top in pyhmmer.hmmer.nhmmer([hmm], [target], cpus=cpus, T=threshold, domT=threshold):
        for hit in top:
            for domain in hit.domains:
                if domain.score < threshold:
                    continue
                aln = domain.alignment
                density = float(domain.score) / (aln.hmm_to - aln.hmm_from + 1)
                step = 1 if aln.target_to >= aln.target_from else -1
                k, t = aln.hmm_from, aln.target_from
                for h, c in zip(aln.hmm_sequence, aln.target_sequence, strict=True):
                    if h != "." and k in by_state and c != "-":
                        for name, orientation in by_state[k]:
                            strand = step * orientation
                            # 1-based forward coordinate of the start-codon base -> 0-based on its strand
                            t_start = t - 1 if strand == 1 else seq_len - t
                            positions.append((float(domain.score), name, strand, t_start, density))
                    if h != ".":
                        k += 1
                    if c != "-":
                        t += step
    return positions


def locus_gap_report(
    dna_seq: str | Seq,
    peptides: list[str] | None = None,
    min_score: float | None = None,
    profiles_path: str | Path | None = None,
    cpus: int = 1,
) -> list[dict[str, Any]]:
    """Peptides whose start-codon column falls in a gap, with the length of that gap.

    A peptide absent from ``locate_mdp_loci`` can mean two different things: the genome has no such
    locus, or the locus is there but the start codon sits in a deletion or an assembly gap. The
    output cannot tell those apart, so this reports the second case explicitly.

    Each entry gives the peptide, the profile that located the locus, the profile column of the
    start codon, and ``gap_columns``: how many consecutive profile columns from there align to no
    target base. Measured in profile columns, not target bases, precisely because the target has
    nothing there to count.

    Not a detection failure: with the codon outside the record there is no frame to read, and
    anchoring on a neighbouring base would fabricate a peptide out of frame. Across 343 primate
    mitogenomes this reports SHLP5 in 57 genomes (16.6%), all with gaps of 8-12 columns.
    """
    seq_str = str(dna_seq).upper()
    profiles = load_locus_profiles(profiles_path)
    alphabet = pyhmmer.easel.Alphabet.dna()
    out: list[dict[str, Any]] = []
    for prof_name, prof in profiles.items():
        wanted = list(prof["peptides"]) if peptides is None else [p for p in peptides if p in prof["peptides"]]
        if not wanted:
            continue
        states = {n: prof["peptides"][n]["state"] for n in wanted}
        threshold = prof["evidence_bitscore"] if min_score is None else min_score
        target = pyhmmer.easel.TextSequence(name=b"target", sequence=seq_str).digitize(alphabet)
        for top in pyhmmer.hmmer.nhmmer([prof["hmm"]], [target], cpus=cpus, T=threshold, domT=threshold):
            for hit in top:
                for domain in hit.domains:
                    if domain.score < threshold:
                        continue
                    aln = domain.alignment
                    column: dict[int, str] = {}
                    k = aln.hmm_from
                    for h, c in zip(aln.hmm_sequence, aln.target_sequence, strict=True):
                        if h != ".":
                            column[k] = c
                            k += 1
                    for name, state in states.items():
                        if column.get(state) != "-":
                            continue
                        run = 0
                        while column.get(state + run) == "-":
                            run += 1
                        out.append(
                            {
                                "name": name,
                                "locus_profile": prof_name,
                                "synteny_score": round(float(domain.score), 4),
                                "start_column": state,
                                "gap_columns": run,
                            }
                        )
    return sorted(out, key=lambda d: (d["name"], d["locus_profile"]))


def locate_mdp_loci(
    dna_seq: str | Seq,
    table: int | str = 2,
    peptides: list[str] | None = None,
    profiles_path: str | Path | None = None,
    min_score: float | None = None,
    max_len: int = 60,
    cpus: int = 1,
) -> list[dict[str, Any]]:
    """Locate mitochondrial-derived peptide loci (Humanin, SHLP1-6) by their position in MT-RNR2.

    Every profile of ``mdp_locus_profiles.json`` (the whole 16S rRNA gene, then the Humanin core) is
    searched with nhmmer on both strands; each hit at or above ``min_score`` (default: the profile's
    evidence threshold, calibrated above the paralogous bacterial/plastid 23S and nuclear 28S rRNAs)
    maps the start codon of each peptide of that profile. The frame is then read as the species' own
    ORF (``_read_locus_frame``): Humanin under ``table``, each SHLP under its primary translation
    table (Cobb et al. 2016). A locus already found by an earlier profile is not repeated. Loci are
    sorted by decreasing ``synteny_score`` (bitscore of the profile that found them, named in
    ``locus_profile``).

    ``locus_density`` is the bitscore per aligned profile column. ``vertebrate_like`` is True when it
    reaches the profile's ``vertebrate_min_density`` (calibrated between non-vertebrate controls and the
    most divergent vertebrate mitogenomes); profiles without that calibration (``Humanin_core``, whose
    density does not separate the groups) give False. See ``mdp_display_name``.
    """
    seq_str = str(dna_seq).upper()
    profiles = load_locus_profiles(str(profiles_path) if profiles_path else None)
    loci: list[dict[str, Any]] = []
    for prof_name, prof in profiles.items():
        peps = {n: p for n, p in prof["peptides"].items() if peptides is None or n in peptides}
        if not peps or len(seq_str) < prof["profile_length"] // 4:
            continue
        threshold = prof["evidence_bitscore"] if min_score is None else min_score
        min_density = prof.get("vertebrate_min_density")
        states = {n: (p["state"], p["orientation"]) for n, p in peps.items()}
        oriented = {1: seq_str}
        for score, name, strand, t_start, density in _profile_positions(
            seq_str, prof["hmm"], states, threshold, cpus
        ):
            if strand not in oriented:
                oriented[strand] = str(Seq(seq_str).reverse_complement())
            pep = peps[name]
            pep_table = pep["primary_translation_table"] or table
            locus = _read_locus_frame(
                oriented[strand],
                t_start,
                strand,
                len(seq_str),
                score,
                pep["reference_orf"],
                pep_table,
                max_len,
                score_pf15040=name == "Humanin",
            )
            if locus is None:
                continue
            if any(
                prev["name"] == name
                and prev["strand"] == strand
                and calculate_overlap(prev["start"], prev["end"], locus["start"], locus["end"]) > 0
                for prev in loci
            ):
                continue
            loci.append(
                {
                    **locus,
                    "name": name,
                    "locus_profile": prof_name,
                    "locus_density": round(density, 4),
                    "vertebrate_like": min_density is not None and density >= min_density,
                }
            )
    loci.sort(key=lambda x: -x["synteny_score"])
    return loci


def mdp_display_name(
    name: str,
    vertebrate_like: bool | None,
    pf15040_supported: bool = False,
    at_locus: bool = True,
) -> str:
    """Reported name of a peptide frame, saying how far the evidence goes.

    Three outcomes, from the two independent routes of evidence (position in an MT-RNR2 locus, and
    the PF15040 Humanin domain at its gathering threshold):

    - ``Humanin`` / ``SHLP1``...: at a vertebrate-like MT-RNR2 locus, or -- for Humanin -- at a locus
      whose frame also carries the PF15040 domain. Both the position and the peptide are supported.
    - ``Humanin_position`` / ``SHLP1_position``...: at an MT-RNR2 locus outside the vertebrate
      pattern. Only the homologous position is marked; no peptide is established there.
    - ``Humanin_domain``: ``at_locus=False``, i.e. a frame carrying the PF15040 domain *outside* any
      located MT-RNR2 locus. The domain is supported, the position is not. In 343 primate mitogenomes
      28 such frames occur under ``scan`` (short frames sharing the ``MITCSLN`` motif), so the name
      keeps them distinct from Humanin proper.

    ``at_locus`` is the position route; ``vertebrate_like`` qualifies it and is only consulted when
    ``at_locus`` is true.

    At a locus, ``vertebrate_like`` must be a measured ``True``/``False``, never ``None``. "Not
    measured" and "measured and below the threshold" are different claims, and treating the first as
    the second is how the human rCRS Humanin once came out named ``Humanin_position`` through the
    removed pairwise route. A caller with no density to report has no basis for either name.
    """
    if not at_locus:
        return f"{name}_domain"
    if vertebrate_like is None:
        raise ValueError(
            f"vertebrate_like is None for {name!r} at a locus: the profile density was not measured, "
            "so neither a plain nor a _position name can be justified."
        )
    if vertebrate_like or (name == "Humanin" and pf15040_supported):
        return name
    return f"{name}_position"


def locate_humanin_loci(
    dna_seq: str | Seq,
    table: int | str = 2,
    profiles_path: str | Path | None = None,
    min_score: float | None = None,
    max_len: int = 60,
    max_hits: int | None = None,
    cpus: int = 1,
) -> list[dict[str, Any]]:
    """Locate Humanin loci by their position in MT-RNR2 and read the species' own frame at each one.

    Uses the multi-clade locus profiles (``locate_mdp_loci``): every mitochondrial and NUMT copy at
    or above ``min_score``. The position only fixes the locus and its homologous reading frame. The
    status is judged on the species' own ORF: the start codon nearest to the mapped human start,
    within ``LOCAL_START_WINDOW`` codons either way, whose ORF covers at least
    ``INTERRUPTION_FRACTION`` of the homologous span (any initiator of a mitochondrial ``table``;
    Met-decoding codons only under other tables). That ORF is ``Complete`` (start through an in-frame
    stop) or ``Incomplete`` (no stop within ``max_len``). Without such an ORF the homologous span read
    from the mapped codon is ``Interrupted`` (stop codon inside it, shown as ``*``) or ``Incomplete``
    (no start codon). ``start_offset_codons`` gives the local start relative to the mapped human one
    (``None`` without a local ORF) and ``identity_to_reference`` the ungapped identity of the
    homologous span to the human peptide; neither enters the status. Loci are sorted by decreasing
    ``synteny_score``.

    A single-sequence rCRS anchor aligned pairwise was once offered here as ``method="pairwise"``.
    It was removed in 1.0.0b1: it could not measure profile density, so it reported
    ``vertebrate_like=None``, and the human rCRS Humanin then came out named ``Humanin_position`` --
    asserting "homologous position only" about the very sequence the threshold was calibrated on.
    Two independent locus finders also meant two things to keep in agreement for no stated benefit.
    """
    loci = locate_mdp_loci(
        dna_seq,
        table=table,
        peptides=["Humanin"],
        profiles_path=profiles_path,
        min_score=min_score,
        max_len=max_len,
        cpus=cpus,
    )
    return loci[:max_hits] if max_hits else loci


def locate_humanin_locus(dna_seq: str | Seq, table: int | str = 2, **kwargs: Any) -> dict[str, Any] | None:
    """Best-scoring Humanin locus of a record (see ``locate_humanin_loci``), or ``None``."""
    loci = locate_humanin_loci(dna_seq, table=table, max_hits=1, **kwargs)
    return loci[0] if loci else None
