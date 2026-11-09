"""Unit tests for humanin_finder.shlp: SHLP references and positional locus annotation."""

from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

from humanin_finder.shlp import annotate_shlps, load_shlp_references

RCRS = Path(__file__).resolve().parents[1] / "examples" / "rcrs_nc_012920_1.fasta"


def test_load_shlp_references():
    data = load_shlp_references()
    assert "shlp_peptides" in data
    peps = data["shlp_peptides"]
    assert len(peps) == 6
    for k in ("SHLP1", "SHLP2", "SHLP3", "SHLP4", "SHLP5", "SHLP6"):
        assert k in peps
        assert len(peps[k]["peptide_sequence_canonical"]) >= 20


def test_annotate_shlps_rcrs_reference():
    full_seq = str(SeqIO.read(str(RCRS), "fasta").seq)
    annotations = {a["name"]: a for a in annotate_shlps(full_seq)}
    assert set(annotations) == {f"SHLP{i}" for i in range(1, 7)}

    # SHLP2 in rCRS: 2091-2168 (1-based) on the reverse strand, plus its stop codon
    shlp2 = annotations["SHLP2"]
    assert shlp2["strand"] == -1
    assert shlp2["end"] == 2168 and shlp2["start"] == 2090 - 3
    assert shlp2["translation_table"] == 1


def test_annotate_shlps_on_reverse_complemented_record():
    full_seq = str(SeqIO.read(str(RCRS), "fasta").seq)
    rc = str(Seq(full_seq).reverse_complement())
    fwd = {a["name"]: a for a in annotate_shlps(full_seq)}
    rev = {a["name"]: a for a in annotate_shlps(rc)}
    for name, a in fwd.items():
        assert rev[name]["strand"] == -a["strand"] and rev[name]["seq"] == a["seq"]
        assert (len(full_seq) - rev[name]["end"], len(full_seq) - rev[name]["start"]) == (a["start"], a["end"])


def test_shipped_table_2_readings_are_the_translation_of_the_reference_orfs():
    """The descriptive Table 2 fields of shlp_references.json must be what the code would translate.

    Until 1.0.0b1 SHLP2 and SHLP4 carried hand-written Table 2 readings that disagreed with the
    actual translation of their reference ORFs (AGA/AGG are stops under the vertebrate mitochondrial code).
    """
    from Bio.Seq import Seq

    from humanin_finder.core import load_locus_profiles

    orfs = load_locus_profiles()["MT-RNR2"]["peptides"]
    for name, meta in load_shlp_references()["shlp_peptides"].items():
        orf = orfs[name]["reference_orf"]
        t1, t2 = str(Seq(orf).translate(table=1)), str(Seq(orf).translate(table=2))
        primary = t1 if meta["primary_translation_table"] == 1 else t2
        assert meta["peptide_sequence_canonical"] == primary, name
        assert meta["mitochondrial_table_2_sequence"] == t2, name
        assert meta["has_mitochondrial_table_2_premature_stop"] == ("*" in t2), name
