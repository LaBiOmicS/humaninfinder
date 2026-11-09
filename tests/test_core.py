"""Unit tests for humanin_finder.core: sORFs, pyhmmer search, translation, and coordinates."""

from humanin_finder.core import (
    calculate_overlap,
    find_sorfs,
    get_translation_table,
    guided_evolutionary_rescue,
    is_degenerate_peptide,
    run_hmm_search_pyhmmer,
    translate_dna,
)

# Reference Humanin DNA coding sequences:
HUMANIN_CDS_T1 = "ATGGCTCCACGAGGGTTCAGCTGTCTCTTACTTTTAACCAGTGAAATTGACCTGCCCGTGAAGAGGCGGGCATAA"
HUMANIN_AA_T1 = "MAPRGFSCLLLLTSEIDLPVKRRA"
HUMANIN_AA_T2 = "MAPRGFSCLLLLTSEIDLPVK"


def test_translation_tables():
    t1 = get_translation_table(1)
    assert "ATG" in t1.start_codons
    assert "TAA" in t1.stop_codons

    t2 = get_translation_table(2)
    assert "ATA" in t2.start_codons
    assert "AGA" in t2.stop_codons
    assert "AGG" in t2.stop_codons


def test_translate_dna():
    assert translate_dna("ATGAAA", table=1) == "MK"
    assert translate_dna("ATGAAAC", table=1) == "MK"
    assert translate_dna("", table=1) == ""
    assert translate_dna("NN", table=1) == ""


def test_find_sorfs_forward():
    flank = "ACGTACGTACGT"
    full_dna = flank + HUMANIN_CDS_T1 + flank
    sorfs = find_sorfs(full_dna, table=1, min_len=10, max_len=30)
    assert len(sorfs) >= 1
    hit = next(s for s in sorfs if s["seq"] == HUMANIN_AA_T1)
    assert hit["strand"] == 1
    assert hit["start"] == len(flank)
    assert hit["end"] == len(flank) + len(HUMANIN_CDS_T1)
    assert hit["status"] == "Complete"


def test_find_sorfs_table2_truncation():
    # In table 2, AGG at codon 22 is a stop codon
    full_dna = "AAAA" + HUMANIN_CDS_T1 + "TTTT"
    sorfs = find_sorfs(full_dna, table=2, min_len=10, max_len=30)
    hit = next((s for s in sorfs if s["seq"] == HUMANIN_AA_T2), None)
    assert hit is not None
    assert hit["length"] == 21
    assert hit.get("table_1_pep") == HUMANIN_AA_T1


def test_overlap():
    assert calculate_overlap(10, 20, 15, 25) == 5
    assert calculate_overlap(10, 20, 20, 30) == 0
    assert calculate_overlap(10, 20, 30, 40) == 0


def test_pyhmmer_in_memory_search():
    dna = "TTTT" + HUMANIN_CDS_T1 + "AAAA"
    hits = run_hmm_search_pyhmmer(dna, table=1)
    assert len(hits) >= 1
    top_hit = hits[0]
    assert top_hit["score"] > 25.0
    assert top_hit["strand"] == 1
    assert "MAPRGFSCLLLLT" in top_hit["matched_pep"]


def test_guided_evolutionary_rescue():
    dna = "TTTT" + HUMANIN_CDS_T1 + "AAAA"
    hits = run_hmm_search_pyhmmer(dna, table=1)
    rescued = guided_evolutionary_rescue(dna, hits, table=1)
    assert len(rescued) >= 1
    r = rescued[0]
    assert "MAPRGFSCLLLL" in r["seq"]
    assert r["rescued"] is True


def test_is_degenerate_peptide():
    assert is_degenerate_peptide("")
    assert is_degenerate_peptide("A" * 15)  # Homopolymeric
    assert is_degenerate_peptide("M***AKL***")  # Excessive stops
    assert not is_degenerate_peptide("MAPRGFSCLLLLTSEIDLPVKRRA")


def test_guided_evolutionary_rescue_pseudogene():
    # Insert internal stop codon (TAG) at codon 11
    pseudo_dna = "TTTT" + HUMANIN_CDS_T1[:33] + "TAG" + HUMANIN_CDS_T1[36:] + "AAAA"
    hits = run_hmm_search_pyhmmer(pseudo_dna, table=1)
    assert len(hits) >= 1
    rescued = guided_evolutionary_rescue(pseudo_dna, hits, table=1)
    assert len(rescued) >= 1
    r = rescued[0]
    assert "*" in r["seq"]
    assert r["status"] == "Interrupted"
