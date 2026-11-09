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
