"""Tests for the coding-evidence controls (``humanin_finder.conservation``)."""

from __future__ import annotations

import json
import random
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq

from humanin_finder import conservation as cons

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"
RCRS = EXAMPLES / "rcrs_nc_012920_1.fasta"


@pytest.fixture(scope="module")
def rcrs() -> str:
    return str(SeqIO.read(str(RCRS), "fasta").seq).upper()


def test_bundled_reference_matches_rcrs_coordinates(rcrs):
    """The bundled reference genes must be exactly the stated rCRS slices."""
    ref = cons.load_reference()
    assert set(ref) == set(cons.GENE_COORDS)
    for gene, (a, b) in cons.GENE_COORDS.items():
        assert ref[gene] == rcrs[a - 1 : b], gene


def test_humanin_window_translates_to_the_human_peptide(rcrs):
    """The window constant must still name the Humanin frame under the vertebrate code."""
    a, b = cons.HUMANIN_WINDOW
    assert str(Seq(rcrs[a - 1 : b]).translate(table=2)) == "MAPRGFSCLLLLTSEIDLPVK"


def test_module_is_independent_of_the_locus_caller():
    """The control must not consult the MT-RNR2 profiles or the Pfam domain.

    This independence is what makes the module a control on the locus interpretation rather than a
    restatement of it, so it is asserted rather than left to convention.
    """
    source = (Path(cons.__file__)).read_text()
    for forbidden in ("MT-RNR2.hmm", "Humanin_core.hmm", "PF15040", "pyhmmer", "locate_mdp_loci", "mdp_locus_profiles"):
        assert forbidden not in source, forbidden
    # and nothing from .core is imported
    assert "from .core" not in source
    assert "import core" not in source


def test_single_genome_cohort_is_identical_to_the_reference(tmp_path, rcrs):
    """A cohort made of the reference alone has no variation and perfect identity."""
    fasta = tmp_path / "one.fasta"
    fasta.write_text(f">rCRS\n{rcrs}\n")
    rec = cons.analyse_cohort(fasta)
    ce = rec["genes"]["MT-RNR2"]["coding_evidence"]
    assert rec["genomes"] == 1
    assert ce["by_frame"]["peptide"]["aa_identity"] == 1.0
    assert ce["by_frame"]["peptide"]["pos3_share"] is None  # no variation at all
    assert ce["by_frame"]["peptide"]["codons"] == 21


def test_positive_control_separates_the_coding_frame(tmp_path, rcrs):
    """On MT-CO1 the test must put the true frame far above the shifted ones.

    Built by mutating third codon positions of MT-CO1, which is what purifying selection on a
    protein leaves behind. If this assertion fails the controls have lost their power and no
    negative result from them means anything.
    """
    a, b = cons.GENE_COORDS["MT-CO1"]
    rng = random.Random(0)
    genomes = []
    for i in range(6):
        seq = list(rcrs)
        for p in range(a - 1 + 2, b, 3):  # third position of each codon
            if rng.random() < 0.3:
                seq[p] = rng.choice([c for c in "ACGT" if c != seq[p]])
        genomes.append(f">g{i}\n{''.join(seq)}\n")
    fasta = tmp_path / "co1.fasta"
    fasta.write_text("".join(genomes))
    rec = cons.analyse_cohort(fasta)
    pc = rec["genes"]["MT-CO1"]["positive_control"]
    assert pc["pos3_share_in_frame_median"] > 0.9
    assert pc["pos3_share_off_frame_median"] < 0.2


def test_rrna_like_variation_does_not_favour_any_frame(tmp_path, rcrs):
    """Variation spread without regard to codons must leave pos3_share near 1/3."""
    a, b = cons.GENE_COORDS["MT-RNR2"]
    rng = random.Random(1)
    genomes = []
    for i in range(8):
        seq = list(rcrs)
        for p in range(a - 1, b):
            if rng.random() < 0.1:
                seq[p] = rng.choice([c for c in "ACGT" if c != seq[p]])
        genomes.append(f">g{i}\n{''.join(seq)}\n")
    fasta = tmp_path / "rrna.fasta"
    fasta.write_text("".join(genomes))
    rec = cons.analyse_cohort(fasta)
    ce = rec["genes"]["MT-RNR2"]["coding_evidence"]
    assert 0.2 < ce["by_frame"]["peptide"]["pos3_share"] < 0.5
    # and the window must not stand out from the gene's own null
    assert 0.05 < ce["pos3_share_rank"] < 0.95


def test_downstream_relaxation_detects_a_real_relaxation(tmp_path, rcrs):
    """Positive control for the relaxation test: if the segment really drifts, it must show.

    Half the genomes get a stop at codon 12 *and* heavy amino-acid change downstream. The test must
    then report higher divergence in the stop-bearing group, which is the signal that peptide-level
    constraint would predict and that the real cohorts do not show.
    """
    a, _b = cons.GENE_COORDS["MT-RNR2"]
    hs = cons.HUMANIN_WINDOW[0] - a  # 0-based offset of the window inside the gene
    rng = random.Random(2)
    genomes = []
    for i in range(10):
        seq = list(rcrs)
        if i < 5:
            stop_at = a - 1 + hs + 3 * 11  # codon 12, 0-based
            seq[stop_at : stop_at + 3] = list("TAA")
            for k in range(12, 21):  # codons 13..21: change first positions
                p = a - 1 + hs + 3 * k
                seq[p] = rng.choice([c for c in "ACGT" if c != seq[p]])
        genomes.append(f">g{i}\n{''.join(seq)}\n")
    fasta = tmp_path / "relax.fasta"
    fasta.write_text("".join(genomes))
    rec = cons.analyse_cohort(fasta)
    rel = rec["genes"]["MT-RNR2"]["downstream_relaxation"]
    assert rel["with_stop"]["n"] == 5
    assert rel["continuous"]["n"] == 5
    assert rel["with_stop"]["mean_aa_div_segment"] > rel["continuous"]["mean_aa_div_segment"]


def test_motif_rank_reports_equal_length_windows(tmp_path, rcrs):
    fasta = tmp_path / "one.fasta"
    fasta.write_text(f">rCRS\n{rcrs}\n")
    rec = cons.analyse_cohort(fasta)
    motifs = rec["genes"]["MT-RNR2"]["coding_evidence"]["motifs"]
    assert motifs["CLLL"]["codons"] == 4
    assert motifs["TSEIDLPVK"]["codons"] == 9


def test_report_and_summary_round_trip(tmp_path, rcrs):
    """A single-genome cohort has no variation, so the summary omits the median-based lines."""
    fasta = tmp_path / "one.fasta"
    fasta.write_text(f">rCRS\n{rcrs}\n")
    rec = cons.analyse_cohort(fasta)
    path = cons.write_report(rec, str(tmp_path / "out"))
    reloaded = json.loads(Path(path).read_text())
    assert reloaded["genomes"] == 1
    text = cons.format_summary(rec)
    assert "pos3_share" in text
    assert "genomes used: 1" in text


def test_summary_reports_the_positive_control_when_there_is_variation(tmp_path, rcrs):
    """With a varying cohort the digest must state the MT-CO1 in-frame/off-frame contrast."""
    a, b = cons.GENE_COORDS["MT-CO1"]
    rng = random.Random(3)
    genomes = []
    for i in range(4):
        seq = list(rcrs)
        for p in range(a - 1 + 2, b, 3):
            if rng.random() < 0.3:
                seq[p] = rng.choice([c for c in "ACGT" if c != seq[p]])
        genomes.append(f">g{i}\n{''.join(seq)}\n")
    fasta = tmp_path / "varying.fasta"
    fasta.write_text("".join(genomes))
    text = cons.format_summary(cons.analyse_cohort(fasta))
    assert "MT-CO1 positive control" in text


def test_empty_cohort_raises(tmp_path):
    fasta = tmp_path / "empty.fasta"
    fasta.write_text("")
    with pytest.raises(ValueError, match="no sequences"):
        cons.analyse_cohort(fasta)


@pytest.mark.parametrize(
    ("label", "body"),
    [
        ("empty record", ">empty\n\n"),
        ("only N", ">n\n" + "N" * 500 + "\n"),
        ("too short", ">short\nACGT\n"),
        ("ambiguity codes", ">iupac\n" + "RYKMSWBDHVN" * 60 + "\n"),
    ],
)
def test_degenerate_input_is_refused_without_a_traceback(tmp_path, label, body):
    """Unusable input must produce a stated result, not a Biopython traceback.

    Biopython raises on an empty target and returns an empty alignment list for one that cannot be
    aligned at all; both used to surface raw (``ValueError: sequence has zero length`` and
    ``IndexError: index out of range``) while ``predict`` handled the same four inputs cleanly.
    """
    fasta = tmp_path / "degenerate.fasta"
    fasta.write_text(body)
    rec = cons.analyse_cohort(fasta)
    host = rec["genes"]["MT-RNR2"]
    assert host["genomes_used"] == 0, label
    assert host["genomes_dropped_low_coverage"] >= 1, label
    text = cons.format_summary(rec)
    assert "no result" in text, label
    assert "covered less than" in text, label


def test_dropped_genomes_are_counted_alongside_usable_ones(tmp_path, rcrs):
    """A cohort mixing usable and unusable records must say how many it set aside."""
    fasta = tmp_path / "mixed.fasta"
    fasta.write_text(f">good1\n{rcrs}\n>good2\n{rcrs}\n>junk\n{'N' * 400}\n")
    rec = cons.analyse_cohort(fasta)
    host = rec["genes"]["MT-RNR2"]
    assert host["genomes_used"] == 2
    assert host["genomes_dropped_low_coverage"] == 1
    assert "genomes dropped: 1" in cons.format_summary(rec)


def test_projection_survives_an_uncountable_number_of_optimal_alignments(monkeypatch):
    """Biopython raises OverflowError on len() when the optimal alignments exceed 2**63.

    This happened on one genome of the 343-primate cohort; the projection only needs the first
    alignment, so it must take it without ever counting them.
    """

    class _Aln:
        aligned = (((0, 4),), ((2, 6),))
        score = 8.0

    class _Uncountable:
        def __len__(self):
            raise OverflowError("number of optimal alignments is larger than 9223372036854775807")

        def __getitem__(self, i):
            raise AssertionError("indexing is not needed for the first alignment")

        def __iter__(self):
            yield _Aln()

    class _Aligner:
        def align(self, ref, target):
            return _Uncountable()

    monkeypatch.setattr(cons, "_aligner", _Aligner())
    assert cons._project("ACGT", "TTACGTTT") == ("ACGT", 8.0)


def test_projection_of_an_unalignable_target_is_all_gaps(monkeypatch):
    class _Aligner:
        def align(self, ref, target):
            return iter(())

    monkeypatch.setattr(cons, "_aligner", _Aligner())
    assert cons._project("ACGT", "ACGT") == ("----", 0.0)


def test_a_zero_statistic_is_ranked_not_dropped():
    """pos3_share of exactly 0 is a measurement (no variation on third positions) and must get a rank.

    A truthiness test once turned it into ``None``, as if the window had not been measured.
    """
    import random as _random

    rng = _random.Random(11)
    ref = "".join(rng.choice("ACGT") for _ in range(300))
    start, width = 60, 63  # window in frame 0 of the gene

    def mutate(seq, positions):
        s = list(seq)
        for p in positions:
            s[p] = next(b for b in "ACGT" if b != s[p])
        return "".join(s)

    outside = [p for p in range(300) if not start <= p < start + width]
    projections = [
        mutate(ref, [start + 3 * k for k in range(0, 21, 2)]),  # first codon positions only
        mutate(ref, rng.sample(outside, 40)),
    ]
    gene = cons.GeneData.build(ref, projections, population=False)
    ce = cons.coding_evidence(gene, (start + 1, start + width), 1)
    assert ce["by_frame"]["peptide"]["pos3_share"] == 0.0
    assert ce["pos3_share_rank"] is not None and 0.0 <= ce["pos3_share_rank"] <= 1.0


def test_cohort_can_be_gzip_compressed(tmp_path, rcrs):
    import gzip as _gzip

    plain, packed = tmp_path / "c.fasta", tmp_path / "c.fasta.gz"
    plain.write_text(f">a\n{rcrs}\n>b\n{rcrs}\n")
    with _gzip.open(packed, "wt") as fh:
        fh.write(plain.read_text())
    assert cons.analyse_cohort(packed) == {**cons.analyse_cohort(plain), "input": str(packed)}
