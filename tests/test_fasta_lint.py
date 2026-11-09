"""Unit tests for FASTA linter in humanin_finder."""

from humanin_finder.fasta_lint import LintOptions, lint_fasta


def test_lint_valid_fasta(tmp_path):
    fasta_file = tmp_path / "valid.fasta"
    fasta_file.write_text(
        ">seq1 Valid record\nATGGCTCCACGAGGGTTCAGCTGTCTCTTACTTTTAACCAGTGAAATTGACCTGCCCGTGAAGAGGCGGGCATAA\n"
    )
    report = lint_fasta(fasta_file, LintOptions(alphabet="dna"))
    assert report["summary"]["valid"] is True
    assert report["summary"]["error_count"] == 0
    assert report["summary"]["record_count"] == 1


def test_lint_invalid_characters(tmp_path):
    fasta_file = tmp_path / "invalid.fasta"
    fasta_file.write_text(">seq1 Corrupt\nATGCZX123\n")
    report = lint_fasta(fasta_file, LintOptions(alphabet="dna"))
    assert report["summary"]["valid"] is False
    assert report["summary"]["error_count"] > 0
    issue_codes = {issue["code"] for issue in report["issues"]}
    assert "alphabet_mismatch" in issue_codes or "invalid_residue" in issue_codes


def test_lint_empty_record(tmp_path):
    fasta_file = tmp_path / "empty.fasta"
    fasta_file.write_text(">empty\n\n")
    report = lint_fasta(fasta_file, LintOptions(alphabet="dna"))
    assert report["summary"]["valid"] is False
    issue_codes = {issue["code"] for issue in report["issues"]}
    assert "empty_sequence" in issue_codes
