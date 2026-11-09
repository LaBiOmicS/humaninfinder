"""Regression tests for defects found in the 2026-10-19 audit (see ARQUITETURA.md §8)."""

import csv
import gzip
import random
from pathlib import Path

import pytest
from Bio import SeqIO
from click.testing import CliRunner

from humanin_finder.cli import filter_non_redundant, main, write_deterministic_outputs
from humanin_finder.core import find_sorfs, get_gathering_threshold, locate_mdp_loci
from humanin_finder.shlp import annotate_shlps, load_shlp_references

RCRS = Path(__file__).resolve().parents[1] / "examples" / "rcrs_nc_012920_1.fasta"
HUMANIN_CDS_T1 = "ATGGCTCCACGAGGGTTCAGCTGTCTCTTACTTTTAACCAGTGAAATTGACCTGCCCGTGAAGAGGCGGGCATAA"


def _random_dna(n: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def test_gathering_threshold_read_from_pfam_profile():
    assert get_gathering_threshold() == 27.0


def test_shlps_of_human_reference_are_complete_in_primary_context():
    seq = str(SeqIO.read(str(RCRS), "fasta").seq)
    refs = load_shlp_references()["shlp_peptides"]
    annotations = {a["name"]: a for a in annotate_shlps(seq)}
    assert set(annotations) == set(refs)
    for name, a in annotations.items():
        # SHLP1/2/4 carry a premature stop under the mitochondrial code, so each SHLP is read in
        # its primary (literature) context, never in the code of the record
        assert a["translation_table"] == refs[name]["primary_translation_table"]
        assert (a["status"], a["seq"], a["identity_to_reference"]) == (
            "Complete",
            refs[name]["peptide_sequence_canonical"],
            1.0,
        )
        assert a["strand"] == refs[name]["strand"] and a["start_offset_codons"] == 0
    # In silico evidence for the locus; the literature tier belongs to the reference peptide
    assert {a["evidence_tier"] for a in annotations.values()} == {"E1"}
    assert annotations["SHLP2"]["reference_evidence_tier"] == "E4"


def test_divergent_shlp_frame_does_not_inherit_the_literature_tier():
    seq = str(SeqIO.read(str(RCRS), "fasta").seq)
    # Non-synonymous change inside SHLP2 (reverse strand, rCRS 2091-2168)
    mutated = seq[:2130] + ("A" if seq[2130] != "A" else "C") + seq[2131:]
    shlp2 = next(a for a in annotate_shlps(mutated) if a["name"] == "SHLP2")
    assert shlp2["seq"] != load_shlp_references()["shlp_peptides"]["SHLP2"]["peptide_sequence_canonical"]
    assert shlp2["evidence_tier"] == "E1" and shlp2["reference_evidence_tier"] is None
    assert shlp2["identity_to_reference"] < 1.0


def test_non_redundant_filter_keeps_antisense_candidates():
    base = {"start_0based": 100, "end_0based": 172, "genetic_code": 2, "hmm_supported": False}
    plus = {**base, "strand": 1, "hmm_bitscore": 40.0}
    minus = {**base, "strand": -1, "hmm_bitscore": 35.0}
    assert len(filter_non_redundant([plus, minus])) == 2
    assert len(filter_non_redundant([plus, {**plus, "hmm_bitscore": 30.0}])) == 1


def test_circular_wrap_coordinates_extend_past_origin():
    rng = random.Random(0)
    body = "".join(rng.choice("CG") for _ in range(300))
    cds = "ATG" + "GCT" * 14 + "TAA"
    circ = cds[30:] + body + cds[:30]
    wraps = [s for s in find_sorfs(circ, table=1, circular=True) if s["circular_wrap"] and s["strand"] == 1]
    assert wraps
    for s in wraps:
        assert 0 <= s["start"] < len(circ) < s["end"]
        assert (circ + circ)[s["start"] : s["end"]] == cds


def test_circular_wrap_minus_strand_maps_to_forward_coordinates():
    rng = random.Random(1)
    body = "".join(rng.choice("CG") for _ in range(300))
    cds = "ATG" + "GCT" * 14 + "TAA"
    circ = cds[30:] + body + cds[:30]
    from Bio.Seq import Seq

    rc = str(Seq(circ).reverse_complement())
    hit = next(s for s in find_sorfs(rc, table=1, circular=True) if s["circular_wrap"] and s["strand"] == -1)
    length = len(rc)
    # Re-extract the forward-strand interval (wrapping) and check it is the reverse complement of the CDS
    fwd = (rc + rc)[hit["start"] : hit["end"]]
    assert str(Seq(fwd).reverse_complement()) == cds
    assert hit["end"] > length


def test_find_sorfs_bounded_search_is_fast_on_long_open_frames():
    import time

    t0 = time.time()
    find_sorfs("ATG" * 20000, table=1, strand_filter=1)
    assert time.time() - t0 < 2.0


def test_bed_splits_origin_crossing_features(tmp_path):
    row = {
        "sample_id": "circ",
        "name": "Humanin",
        "status": "Canonical",
        "seq": "M" * 16,
        "strand": 1,
        "start_0based": 318,
        "end_0based": 369,
        "start_1based": 319,
        "end_1based": 369,
        "seq_length": 360,
        "seq_circular": True,
    }
    write_deterministic_outputs([row], str(tmp_path / "c"))
    bed = [line.split("\t") for line in (tmp_path / "c_results.bed").read_text().splitlines() if line[0] != "#"]
    assert [(b[1], b[2]) for b in bed] == [("318", "360"), ("0", "9")]
    gff = (tmp_path / "c_results.gff3").read_text()
    assert "Is_circular=true" in gff and "##sequence-region circ 1 360" in gff


def test_classifier_options_are_gone(tmp_path):
    for opt in (["-t", "0.5"], ["--classifier"], ["--device", "cpu"], ["--batch-size", "8"]):
        result = CliRunner().invoke(main, ["predict", "-i", str(RCRS), "-o", str(tmp_path / "x"), *opt])
        assert result.exit_code != 0 and "No such option" in result.output


def test_no_hmm_on_large_sequence_warns(tmp_path):
    fa = tmp_path / "big.fasta"
    fa.write_text(">big\n" + _random_dna(60_000, 5) + "\n")
    result = CliRunner().invoke(
        main, ["scan", "-i", str(fa), "-o", str(tmp_path / "b"), "--no-hmm", "--no-synteny", "-g", "1"]
    )
    assert result.exit_code == 0, result.output
    assert "enable --synteny or --hmm" in result.output


def test_gzip_input_runs_end_to_end(tmp_path):
    gz = tmp_path / "rcrs.fasta.gz"
    with gzip.open(gz, "wt") as handle:
        handle.write(RCRS.read_text())
    result = CliRunner().invoke(main, ["predict", "-i", str(gz), "-o", str(tmp_path / "gz")])
    assert result.exit_code == 0, result.output
    assert "2633-2698" in result.output


def test_scan_reports_every_numt_like_locus(tmp_path):
    seq = _random_dna(40_000, 3) + HUMANIN_CDS_T1 + _random_dna(60_000, 4) + HUMANIN_CDS_T1 + _random_dna(40_000, 6)
    fa = tmp_path / "numt.fasta"
    fa.write_text(">chr_syn\n" + seq + "\n")
    out = tmp_path / "numt"
    result = CliRunner().invoke(main, ["scan", "-i", str(fa), "-o", str(out), "-g", "1", "--no-shlp"])
    assert result.exit_code == 0, result.output
    rows = _read_csv(f"{out}_results.csv")
    hmm_loci = [r for r in rows if r["hmm_supported"] == "True"]
    # Each copy is identified by its stop codon: an in-frame upstream ATG in the random flank may extend the ORF
    n = len(HUMANIN_CDS_T1)
    assert {40_000 + n, 40_000 + n + 60_000 + n} <= {int(r["end_1based"]) for r in hmm_loci}


def test_outputs_are_honest_about_evidence(tmp_path):
    out = tmp_path / "rcrs"
    result = CliRunner().invoke(main, ["predict", "-i", str(RCRS), "-o", str(out), "--all-candidates"])
    assert result.exit_code == 0, result.output
    rows = _read_csv(f"{out}_results.csv")
    assert "prob_calibrated" not in rows[0]
    shlp = [r for r in rows if r["name"].startswith("SHLP")]
    assert {r["status"] for r in shlp} == {"Complete"} and len(shlp) == 6
    assert all(r["hmm_evalue"] == "" for r in shlp)
    assert all(float(r["locus_density"]) > 0.8 for r in shlp)
    assert {r["evidence_tier"] for r in rows} == {"E1"}
    assert not {"classifier_score", "in_training_domain"} & set(rows[0])
    humanin = [r for r in rows if r["name"] == "Humanin"]
    # Humanin is named by homology only (MT-RNR2 position and/or PF15040)
    assert all("Synteny" in r["evidence_source"] or "HMM" in r["evidence_source"] for r in humanin)
    assert {int(r["start_1based"]) for r in humanin} == {2633}
    assert not any(r["name"].startswith("HLP") or "_position" in r["name"] for r in rows)
    assert {r["status"] for r in rows} <= {"Complete", "Incomplete", "Interrupted"}
    gff = Path(f"{out}_results.gff3").read_text()
    assert "\tORF\t" in gff and "\tpolypeptide\t" not in gff


def test_alternative_start_codon_is_decoded_as_initiator_met():
    # GTG is a vertebrate mitochondrial start codon (Table 2) that reads as Val internally
    orf = "GTG" + "GCT" * 12 + "TAA"
    hit = next(s for s in find_sorfs("CC" + orf + "CC", table=2, allow_alternative_starts=True) if s["start"] == 2)
    assert hit["seq"].startswith("M") and hit["start_codon"] == "GTG"


def test_table_1_reading_requires_a_standard_start_codon():
    ata_orf = "ATA" + "GCT" * 12 + "TAA"  # ATA initiates under Table 2 but not under Table 1
    atg_orf = "ATG" + "GCT" * 12 + "TAA"
    ata = next(s for s in find_sorfs("CC" + ata_orf + "CC", table=2) if s["start"] == 2)
    atg = next(s for s in find_sorfs("CC" + atg_orf + "CC", table=2) if s["start"] == 2)
    assert "table_1_pep" not in ata
    assert atg["table_1_pep"].startswith("M")


def test_interrupted_frames_are_reported_by_hmm_support(tmp_path):
    pseudo = HUMANIN_CDS_T1[:33] + "TAG" + HUMANIN_CDS_T1[36:]
    fa = tmp_path / "pseudo.fasta"
    fa.write_text(">pseudo\n" + _random_dna(300, 8) + pseudo + _random_dna(300, 9) + "\n")
    out = tmp_path / "pseudo"
    result = CliRunner().invoke(
        main, ["predict", "-i", str(fa), "-o", str(out), "-g", "1", "--all-candidates"]
    )
    assert result.exit_code == 0, result.output
    rows = _read_csv(f"{out}_results.csv")
    interrupted = [r for r in rows if r["status"] == "Interrupted"]
    assert len(interrupted) >= 1
    assert all(r["hmm_supported"] == "True" for r in interrupted)
    assert {r["evidence_source"] for r in interrupted} == {"HMM"}


def test_unwrapped_records_longer_than_one_megabase_pass_pipeline_lint(tmp_path):
    from humanin_finder.cli import PIPELINE_LINT_OPTIONS
    from humanin_finder.fasta_lint import lint_fasta

    fa = tmp_path / "unwrapped.fasta"
    fa.write_text(">chr\n" + "ACGT" * 300_000 + "\n")  # 1.2 Mb on a single line
    assert lint_fasta(fa, PIPELINE_LINT_OPTIONS)["summary"]["valid"]


def _read_csv(path: str) -> list[dict[str, str]]:
    """Rows of an output CSV as plain strings, skipping the provenance comment lines."""
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(line for line in fh if not line.startswith("#")))


def _rcrs() -> str:
    return str(SeqIO.read(str(RCRS), "fasta").seq).upper()


def test_synteny_locus_recovers_human_humanin_in_both_codes():
    from humanin_finder.core import locate_humanin_locus

    seq = _rcrs()
    t2, t1 = locate_humanin_locus(seq, table=2), locate_humanin_locus(seq, table=1)
    assert (t2["seq"], t2["start"], t2["end"], t2["status"]) == ("MAPRGFSCLLLLTSEIDLPVK", 2632, 2698, "Complete")
    assert (t1["seq"], t1["start"], t1["end"], t1["status"]) == ("MAPRGFSCLLLLTSEIDLPVKRRA", 2632, 2707, "Complete")


def test_synteny_locus_maps_reverse_strand_to_forward_coordinates():
    from Bio.Seq import Seq

    from humanin_finder.core import locate_humanin_locus

    seq = _rcrs()
    rc = str(Seq(seq).reverse_complement())
    hit = locate_humanin_locus(rc, table=2)
    assert hit["strand"] == -1 and (len(seq) - hit["end"], len(seq) - hit["start"]) == (2632, 2698)


def test_synteny_locus_absent_in_shuffled_genome():
    from humanin_finder.core import locate_humanin_locus

    shuffled = list(_rcrs())
    random.Random(3).shuffle(shuffled)
    assert locate_humanin_locus("".join(shuffled), table=2) is None


def test_synteny_locus_reports_interrupted_and_incomplete_frames():
    from humanin_finder.core import locate_humanin_locus

    seq = _rcrs()
    interrupted = seq[:2638] + "TAA" + seq[2641:]  # stop at codon 3, as in teleost loci
    no_start = seq[:2632] + "ACG" + seq[2635:]  # mapped codon no longer initiates
    hit_i = locate_humanin_locus(interrupted, table=2)
    hit_n = locate_humanin_locus(no_start, table=2)
    assert hit_i["status"] == "Interrupted" and hit_i["seq"].startswith("MA*")
    assert hit_n["status"] == "Incomplete" and hit_n["seq"].startswith("T")


def test_synteny_status_follows_the_local_orf_not_the_human_start():
    from humanin_finder.core import locate_humanin_locus

    seq = _rcrs()
    # Mapped codon no longer initiates, but the species' own start lies one codon downstream
    shifted = seq[:2632] + "ACGATG" + seq[2638:]
    hit = locate_humanin_locus(shifted, table=2)
    assert (hit["status"], hit["seq"], hit["start"], hit["start_offset_codons"]) == (
        "Complete",
        "MPRGFSCLLLLTSEIDLPVK",
        2635,
        1,
    )
    human = locate_humanin_locus(seq, table=2)
    assert human["start_offset_codons"] == 0 and human["identity_to_reference"] == 1.0
    # A start too far downstream to cover the homologous span leaves the frame interrupted
    interrupted = seq[:2638] + "TAA" + seq[2641:]
    assert locate_humanin_locus(interrupted, table=2)["start_offset_codons"] is None


def test_the_locus_route_is_the_only_one_and_measures_density():
    """One locus finder, and it always measures the density the naming depends on.

    A second route (a single rCRS anchor aligned pairwise) was removed in 1.0.0b1: it could not
    measure profile density, so it reported ``vertebrate_like=None`` and the human rCRS Humanin came
    out named ``Humanin_position``. The signature no longer accepts ``method``, and density is always
    a measured boolean.
    """
    import pytest

    from humanin_finder.core import locate_humanin_locus, mdp_display_name

    seq = _rcrs()
    hit = locate_humanin_locus(seq, table=2)
    assert hit["locus_profile"] == "MT-RNR2"
    assert hit["synteny_score"] >= 40
    assert hit["vertebrate_like"] is True, "the reference must be measured, not assumed"
    assert isinstance(hit["locus_density"], float)
    assert mdp_display_name(hit["name"], hit["vertebrate_like"]) == "Humanin"

    with pytest.raises(TypeError):
        locate_humanin_locus(seq, table=2, method="pairwise")
    # and an unmeasured density is refused rather than silently named _position
    with pytest.raises(ValueError, match="was not measured"):
        mdp_display_name("Humanin", None)


def test_locus_profile_reports_every_copy_in_a_record():
    from Bio.Seq import Seq

    from humanin_finder.core import locate_humanin_loci

    seq = _rcrs()
    segment = seq[2300:3100]
    spacer = "".join(random.Random(5).choice("ACGT") for _ in range(5000))
    record = spacer + segment + spacer + str(Seq(segment).reverse_complement()) + spacer
    loci = locate_humanin_loci(record, table=1)
    assert sorted((locus["strand"], locus["seq"]) for locus in loci) == [
        (-1, "MAPRGFSCLLLLTSEIDLPVKRRA"),
        (1, "MAPRGFSCLLLLTSEIDLPVKRRA"),
    ]


def test_position_only_names_outside_vertebrate_like_loci():
    from humanin_finder.core import mdp_display_name

    assert mdp_display_name("Humanin", True) == "Humanin"
    assert mdp_display_name("Humanin", False) == "Humanin_position"
    # A PF15040 domain at GA names the frame Humanin whatever the locus
    assert mdp_display_name("Humanin", False, pf15040_supported=True) == "Humanin"
    assert mdp_display_name("SHLP3", False, pf15040_supported=True) == "SHLP3_position"
    # An unmeasured density used to be named _position here, which asserts "homologous position only"
    # without having looked. It is now refused; see test_the_locus_route_is_the_only_one_and_...
    with pytest.raises(ValueError, match="was not measured"):
        mdp_display_name("Humanin", None)


def test_humanin_core_loci_are_not_vertebrate_like():
    seq = str(SeqIO.read(str(RCRS), "fasta").seq)
    # The 211-nt core alone: found only by Humanin_core, whose density is not calibrated
    core = seq[2572:2783]
    loci = locate_mdp_loci(core, table=2, peptides=["Humanin"])
    assert loci and all(lc["locus_profile"] == "Humanin_core" and lc["vertebrate_like"] is False for lc in loci)
    full = locate_mdp_loci(seq, table=2, peptides=["Humanin"])
    assert full[0]["locus_profile"] == "MT-RNR2" and full[0]["vertebrate_like"] is True


def test_auto_genetic_code_by_record_type():
    from humanin_finder.cli import auto_genetic_code, parse_genetic_codes

    assert parse_genetic_codes("auto") == []
    assert auto_genetic_code("NC_012920.1 Homo sapiens mitochondrion, complete genome", 16569) == 2
    assert auto_genetic_code("contig_7", 30_000) == 2
    assert auto_genetic_code("NC_000017.11 Homo sapiens chromosome 17", 83_257_441) == 1


def test_local_start_outside_mitochondria_requires_met():
    from humanin_finder.core import locate_humanin_locus

    seq = _rcrs()
    # Mapped ATG removed; a CTG (near-cognate, Table 1 start codon list) lies one codon downstream
    ctg = seq[:2632] + "ACGCTG" + seq[2638:]
    assert locate_humanin_locus(ctg, table=1)["start_offset_codons"] is None
    atg = seq[:2632] + "ACGATG" + seq[2638:]
    assert locate_humanin_locus(atg, table=1)["start_offset_codons"] == 1


def test_gapped_start_codon_is_reported_as_a_gap_not_as_a_peptide():
    """A start codon inside a deletion yields no frame, and that fact must be reportable.

    SHLP5 is absent from 57 of 343 primate mitogenomes because its start-codon column aligns to a
    gap of 8-12 profile columns: the codon is not in the record. Anchoring on a neighbouring base
    instead fabricates an out-of-frame peptide (identity to the human peptide fell to 0.04, with
    internal stops), so ``locate_mdp_loci`` stays silent and ``locus_gap_report`` says why.

    Exercised here on Humanin, whose start codon the reference places at rCRS 2633, so that the case
    is reproducible from the bundled rCRS alone.
    """
    from humanin_finder.core import locate_mdp_loci, locus_gap_report

    seq = _rcrs()
    assert locate_mdp_loci(seq, table=2)
    assert locus_gap_report(seq) == [], "the reference itself must report no gaps"

    deleted = seq[:2630] + seq[2645:]  # removes the Humanin start codon
    names = {lc["name"] for lc in locate_mdp_loci(deleted, table=2)}
    report = locus_gap_report(deleted)
    assert "Humanin" not in names, "an unanchored start codon must not yield a peptide"
    assert report, "and the gap must be reported instead of passing in silence"
    assert {d["name"] for d in report} == {"Humanin"}
    # both profiles that carry Humanin report it, each with its own column
    assert {d["locus_profile"] for d in report} == {"MT-RNR2", "Humanin_core"}
    assert all(d["gap_columns"] >= 1 for d in report)
    # the SHLP loci, whose start codons were not touched, are still found
    assert {"SHLP1", "SHLP2", "SHLP3", "SHLP4", "SHLP5", "SHLP6"} <= names



def test_domain_without_locus_is_named_humanin_domain():
    """PF15040 at the gathering threshold outside any MT-RNR2 locus is not plain ``Humanin``.

    Symmetric to ``Humanin_position`` (position without peptide): this is peptide without position.
    """
    from humanin_finder.core import mdp_display_name

    assert mdp_display_name("Humanin", None, pf15040_supported=True, at_locus=False) == "Humanin_domain"
    assert mdp_display_name("Humanin", True, pf15040_supported=True, at_locus=False) == "Humanin_domain"
    # at a locus the previous behaviour is unchanged
    assert mdp_display_name("Humanin", False, pf15040_supported=True, at_locus=True) == "Humanin"
    assert mdp_display_name("Humanin", True, pf15040_supported=False, at_locus=True) == "Humanin"
