"""Unit and integration tests for humanin-finder CLI."""

from pathlib import Path

from click.testing import CliRunner

from humanin_finder import __version__
from humanin_finder.cli import main, write_deterministic_outputs


def test_cli_version():
    runner = CliRunner()
    result = runner.invoke(main, ["--version"])
    assert result.exit_code == 0
    assert __version__ in result.output


def test_cli_lint():
    runner = CliRunner()
    fasta_path = Path(__file__).resolve().parents[1] / "examples" / "rcrs_nc_012920_1.fasta"
    result = runner.invoke(main, ["lint", str(fasta_path)])
    assert result.exit_code == 0
    assert "100% valid" in result.output


def test_write_deterministic_outputs(tmp_path):
    out_prefix = str(tmp_path / "test_out")
    dummy_results = [
        {
            "sample_id": "seq1",
            "name": "Humanin",
            "status": "Complete",
            "seq": "MAPRGFSCLLLLTSEIDLPVKRRA",
            "length": 24,
            "strand": 1,
            "frame": 0,
            "start_0based": 2632,
            "end_0based": 2704,
            "start_1based": 2633,
            "end_1based": 2704,
            "hmm_supported": True,
            "evidence_source": "HMM",
            "hmm_bitscore": 58.2,
            "hmm_evalue": 1.2e-15,
            "genetic_code": 2,
            "locus_mode": "Scan",
            "evidence_tier": "E1",
            "table_1_pep": "MAPRGFSCLLLLTSEIDLPVKRRA",
        }
    ]

    write_deterministic_outputs(dummy_results, out_prefix)

    csv_file = Path(f"{out_prefix}_results.csv")
    fasta_file = Path(f"{out_prefix}_results.fasta")
    gff3_file = Path(f"{out_prefix}_results.gff3")
    bed_file = Path(f"{out_prefix}_results.bed")

    assert csv_file.exists()
    assert fasta_file.exists()
    assert gff3_file.exists()
    assert bed_file.exists()

    with open(gff3_file) as f:
        gff_lines = f.readlines()
    assert any("HumaninFinder\tORF" in line for line in gff_lines)

    with open(bed_file) as f:
        bed_lines = f.readlines()
    assert any("2632\t2704" in line for line in bed_lines)
    # BED score: 1000 for named loci, 500 for loci that only mark a homologous position
    assert any(line.split("\t")[4] == "1000" for line in bed_lines if not line.startswith("#"))
    write_deterministic_outputs([{**dummy_results[0], "name": "Humanin_position"}], out_prefix + "_pos")
    bed_pos = Path(f"{out_prefix}_pos_results.bed").read_text().splitlines()
    assert [line.split("\t")[4] for line in bed_pos if not line.startswith("#")] == ["500"]


def test_parse_genetic_codes():
    from humanin_finder.cli import parse_genetic_codes

    assert parse_genetic_codes("2") == [2]
    assert parse_genetic_codes("2,1") == [2, 1]
    assert parse_genetic_codes(" 2 , 1 , 5 ") == [2, 1, 5]
    assert parse_genetic_codes(2) == [2]


def test_cli_multi_table(tmp_path):
    runner = CliRunner()
    fasta_path = Path(__file__).resolve().parents[1] / "examples" / "rcrs_nc_012920_1.fasta"
    out_prefix = str(tmp_path / "multi_test")
    result = runner.invoke(main, ["predict", "-i", str(fasta_path), "-o", out_prefix, "-g", "2,1"])
    assert result.exit_code == 0
    assert "Humanin_T2" in result.output
    assert "Humanin_T1" in result.output
    assert "SHLP1" in result.output


def test_cli_large_genome_locus_zoom(tmp_path):
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    from tests.test_core import HUMANIN_CDS_T1

    runner = CliRunner()
    # 200 kb chromosome with Humanin at 120,000
    seq = ("A" * 60000 + "C" * 60000) + HUMANIN_CDS_T1 + ("G" * 80000)
    rec = SeqRecord(Seq(seq), id="chr_syn_200k", description="Synthetic 200kb chromosome")
    test_fa = tmp_path / "syn_200k.fasta"
    SeqIO.write([rec], test_fa, "fasta")

    out_prefix = str(tmp_path / "zoom_out")
    result = runner.invoke(main, ["predict", "-i", str(test_fa), "-o", out_prefix, "-g", "1"])
    assert result.exit_code == 0
    assert "chr_syn_200k" in result.output
    assert "120001-120075" in result.output
    # No 16S rRNA context: the PF15040 domain holds but the position does not, so the frame is
    # named Humanin_domain rather than Humanin (see core.mdp_display_name)
    assert "Humanin_domain" in result.output


def test_write_csv_keeps_integer_columns_integer(tmp_path):
    """An integer column with empty cells stays integer (pandas wrote ``-1.0``; fixed in 1.0.0b1)."""
    from humanin_finder.cli import write_csv

    rows = [
        {"start_offset_codons": -1, "score": 30, "ok": True},
        {"start_offset_codons": None, "score": 12.5, "ok": False},
        {"start_offset_codons": 0, "score": None, "ok": None},
    ]
    out = tmp_path / "t.csv"
    write_csv(rows, ["start_offset_codons", "score", "ok"], str(out))
    assert out.read_text() == "start_offset_codons,score,ok\n-1,30.0,True\n,12.5,False\n0,,\n"
