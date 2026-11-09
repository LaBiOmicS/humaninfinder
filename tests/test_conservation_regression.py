"""Numerical regression tests freezing the conservation numbers that the paper reports.

The controls in ``humanin_finder.conservation`` project each genome onto the reference genes with
``Bio.Align.PairwiseAligner`` and take the first optimal alignment: when several score equally, which
one comes first is an implementation detail of Biopython. Everything downstream -- the codon-position
share, the window's rank inside the 16S null, the motif ranks -- is computed from that projection.

Property tests (``test_conservation.py``) check that the controls still have power: the MT-CO1
positive control must separate the true frame, unstructured rRNA variation must leave pos3_share near
1/3. They would all still pass if the projection shifted and every published figure moved with it.
These tests close that gap by pinning the figures themselves.

If one of them fails, the numbers in the manuscript no longer match what the code produces. That is
not a test to relax: find what changed (a dependency, the aligner parameters, the window constants),
decide whether the new numbers are the correct ones, and only then re-freeze them here, updating the
manuscript and ``tests/data/primates343_conservation_reference.json`` in the same commit.
"""

from __future__ import annotations

import json
import os
import random
from pathlib import Path

import pytest
from Bio import SeqIO

from humanin_finder import conservation as cons

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"
RCRS = EXAMPLES / "rcrs_nc_012920_1.fasta"
REFERENCE = Path(__file__).resolve().parent / "data" / "primates343_conservation_reference.json"

# Tight enough to catch any change of algorithm, loose enough for platform float noise.
REL = 1e-9

# Values produced by HumaninFinder 1.0.0b1 with biopython 1.88 on the cohort built by
# ``synthetic_cohort`` below. Frozen deliberately; see the module docstring before changing them.
FROZEN = {
    "genomes": 12,
    "peptide": {
        "codons": 21,
        "aa_identity": 0.8571428571428571,
        "aa_invariant_codons": 3,
        "pos3_share": 0.4117647058823529,
        "var_per_site": 0.06746031746031746,
        "pos_var": [0.8333, 1.6667, 1.75],
    },
    "shifted_1": {
        "codons": 20,
        "aa_identity": 0.8208333333333333,
        "aa_invariant_codons": 1,
        "pos3_share": 0.196078431372549,
    },
    "shifted_2": {
        "codons": 20,
        "aa_identity": 0.8625,
        "aa_invariant_codons": 6,
        "pos3_share": 0.392156862745098,
    },
    "pos3_share_rank": 0.8984635938543755,
    "var_per_site_rank": 0.41282565130260523,
    "motif_CLLL": {
        "codons": 4,
        "aa_identity": 0.8958333333333333,
        "aa_identity_rank": 0.83984375,
        "invariant_fraction": 0.0,
        "invariant_fraction_rank": 0.580078125,
        "null_aa_identity_median": 0.8333333333333334,
        "null_windows": 512,
        "fully_invariant_windows": 0,
    },
    "positive_control": {"in_frame": 1.0, "off_frame": 0.0},
}


@pytest.fixture(scope="module")
def rcrs() -> str:
    return str(SeqIO.read(str(RCRS), "fasta").seq).upper()


def synthetic_cohort(rcrs: str) -> str:
    """Twelve genomes from one fixed seed: rRNA drift, CO1 third positions, and two indels.

    The indels matter: they are what forces the aligner to choose, which is exactly the step whose
    stability these tests exist to guard. Must stay byte-for-byte reproducible, so neither the seed
    nor the order of the random draws may change.
    """
    rng = random.Random(20261109)
    a, b = cons.GENE_COORDS["MT-RNR2"]
    c, d = cons.GENE_COORDS["MT-CO1"]
    out = [f">ref\n{rcrs}\n"]
    for i in range(11):
        s = list(rcrs)
        for p in range(a - 1, b):
            if rng.random() < 0.08:
                s[p] = rng.choice([x for x in "ACGT" if x != s[p]])
        for p in range(c - 1 + 2, d, 3):
            if rng.random() < 0.25:
                s[p] = rng.choice([x for x in "ACGT" if x != s[p]])
        g = "".join(s)
        if i % 4 == 1:
            g = g[: a + 200] + g[a + 203 :]
        if i % 4 == 2:
            g = g[: a + 400] + "ACGTA" + g[a + 400 :]
        out.append(f">g{i}\n{g}\n")
    return "".join(out)


def test_the_synthetic_cohort_is_itself_reproducible(rcrs):
    """If the cohort drifts, the frozen numbers below mean nothing."""
    assert synthetic_cohort(rcrs) == synthetic_cohort(rcrs)
    assert synthetic_cohort(rcrs).count(">") == 12


@pytest.fixture(scope="module")
def analysed(tmp_path_factory, rcrs):
    path = tmp_path_factory.mktemp("golden") / "cohort.fasta"
    path.write_text(synthetic_cohort(rcrs))
    return cons.analyse_cohort(path)


def test_frozen_codon_position_numbers(analysed):
    ce = analysed["genes"]["MT-RNR2"]["coding_evidence"]
    assert analysed["genomes"] == FROZEN["genomes"]
    for frame in ("peptide", "shifted_1", "shifted_2"):
        got, want = ce["by_frame"][frame], FROZEN[frame]
        for key, value in want.items():
            if isinstance(value, list):
                assert got[key] == pytest.approx(value, rel=REL), f"{frame}.{key}"
            elif isinstance(value, float):
                assert got[key] == pytest.approx(value, rel=REL), f"{frame}.{key}"
            else:
                assert got[key] == value, f"{frame}.{key}"


def test_frozen_null_ranks(analysed):
    ce = analysed["genes"]["MT-RNR2"]["coding_evidence"]
    assert ce["pos3_share_rank"] == pytest.approx(FROZEN["pos3_share_rank"], rel=REL)
    assert ce["var_per_site_rank"] == pytest.approx(FROZEN["var_per_site_rank"], rel=REL)


def test_frozen_motif_numbers(analysed):
    got, want = analysed["genes"]["MT-RNR2"]["coding_evidence"]["motifs"]["CLLL"], FROZEN["motif_CLLL"]
    for key, value in want.items():
        if isinstance(value, float):
            assert got[key] == pytest.approx(value, rel=REL), key
        else:
            assert got[key] == value, key


def test_frozen_positive_control(analysed):
    pc = analysed["genes"]["MT-CO1"]["positive_control"]
    assert pc["pos3_share_in_frame_median"] == pytest.approx(FROZEN["positive_control"]["in_frame"], rel=REL)
    assert pc["pos3_share_off_frame_median"] == pytest.approx(FROZEN["positive_control"]["off_frame"], rel=REL)


def test_the_published_primate_cohort_still_gives_the_published_numbers():
    """The real check, when the cohort is at hand: 343 primate mitogenomes, as reported.

    The FASTA is far too large to ship, so point ``HUMANINFINDER_PRIMATES343`` at it to run this.
    Without it the test skips -- which is why the synthetic cohort above is not optional.
    """
    path = os.environ.get("HUMANINFINDER_PRIMATES343")
    if not path or not Path(path).exists():
        pytest.skip("set HUMANINFINDER_PRIMATES343 to the 343-primate FASTA to run this")
    want = json.loads(REFERENCE.read_text())
    rec = cons.analyse_cohort(Path(path), workers=max(1, (os.cpu_count() or 2) // 2))
    ce = rec["genes"]["MT-RNR2"]["coding_evidence"]

    assert rec["genomes"] == want["genomes"]
    assert rec["mode"] == want["mode"]
    for frame, expected in want["by_frame"].items():
        for key, value in expected.items():
            got = ce["by_frame"][frame][key]
            if isinstance(value, (int, float)) and not isinstance(value, bool):
                assert got == pytest.approx(value, rel=REL), f"{frame}.{key}"
            else:
                assert got == value, f"{frame}.{key}"
    assert ce["pos3_share_rank"] == pytest.approx(want["pos3_share_rank"], rel=REL)
    for motif, expected in want["motifs"].items():
        for key, value in expected.items():
            got = ce["motifs"][motif][key]
            if isinstance(value, (int, float)) and not isinstance(value, bool):
                assert got == pytest.approx(value, rel=REL), f"{motif}.{key}"
            else:
                assert got == value, f"{motif}.{key}"
    pc = rec["genes"]["MT-CO1"]["positive_control"]
    assert pc["pos3_share_in_frame_median"] == pytest.approx(want["positive_control"]["in_frame"], rel=REL)
    assert pc["pos3_share_off_frame_median"] == pytest.approx(want["positive_control"]["off_frame"], rel=REL)
    rel_got = rec["genes"]["MT-RNR2"]["downstream_relaxation"]
    assert rel_got["with_stop"]["n"] == want["downstream_relaxation"]["with_stop_n"]
    assert rel_got["continuous"]["n"] == want["downstream_relaxation"]["continuous_n"]
