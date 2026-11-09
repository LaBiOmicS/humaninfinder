"""Strict, non-mutating FASTA input validation and reporting for HumaninFinder."""

from __future__ import annotations

import gzip
import hashlib
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Literal

from Bio import SeqIO

Alphabet = Literal["auto", "dna", "rna", "protein"]
Severity = Literal["error", "warning"]

KNOWN_FASTA_SUFFIXES = frozenset({".fa", ".fasta", ".fas", ".fsa", ".fna", ".ffn", ".faa", ".seq"})
DNA_IUPAC = frozenset("ACGTRYSWKMBDHVN")
RNA_IUPAC = frozenset("ACGURYSWKMBDHVN")
PROTEIN_IUPAC = frozenset("ACDEFGHIKLMNPQRSTVWYBXZJUO*")
NUCLEOTIDE_AMBIGUITY = frozenset("RYSWKMBDHVN")
PROTEIN_AMBIGUITY = frozenset("BXZJUO*")
DEFAULT_MAX_INPUT_BYTES = 2 * 1024 * 1024 * 1024
DEFAULT_MAX_RECORDS = 1_000_000
DEFAULT_MAX_LINE_LENGTH = 1_000_000
DEFAULT_MAX_SEQUENCE_LENGTH = 50_000_000
DEFAULT_MAX_TOTAL_SEQUENCE_LENGTH = 2_000_000_000


@dataclass(frozen=True)
class Issue:
    """One stable machine-readable lint finding."""

    severity: Severity
    code: str
    message: str
    line: int | None = None
    record_id: str | None = None


@dataclass(frozen=True)
class RecordSummary:
    """Summary of one record without changing its source representation."""

    record_id: str
    description: str
    sequence_length: int
    inferred_alphabet: str
    sha256: str


@dataclass(frozen=True)
class PreflightRecord:
    """Source-sequence metadata used to verify the Biopython parse."""

    record_id: str
    description: str
    header_line: int
    sequence_length: int
    sequence_sha256: str


@dataclass(frozen=True)
class LintOptions:
    """Explicit input policy for one lint run."""

    alphabet: Alphabet = "auto"
    disallow_ambiguity: bool = False
    max_input_bytes: int = DEFAULT_MAX_INPUT_BYTES
    max_records: int = DEFAULT_MAX_RECORDS
    max_line_length: int = DEFAULT_MAX_LINE_LENGTH
    max_sequence_length: int = DEFAULT_MAX_SEQUENCE_LENGTH
    max_total_sequence_length: int = DEFAULT_MAX_TOTAL_SEQUENCE_LENGTH

    def __post_init__(self) -> None:
        if self.alphabet not in {"auto", "dna", "rna", "protein"}:
            raise ValueError(f"Unsupported alphabet policy: {self.alphabet!r}.")
        for field_name in (
            "max_input_bytes",
            "max_records",
            "max_line_length",
            "max_sequence_length",
            "max_total_sequence_length",
        ):
            value = getattr(self, field_name)
            if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
                raise ValueError(f"{field_name} must be a positive integer.")


class ResourceLimitExceeded(Exception):
    """A configured input limit prevented further processing."""

    def __init__(self, code: str, message: str) -> None:
        super().__init__(message)
        self.code = code


def open_fasta_text(path: str | Path):
    """Open a plain or gzip-compressed FASTA file in text mode (gzip detected by its magic bytes)."""
    with open(path, "rb") as handle:
        is_gzip = handle.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if is_gzip else open(path)


def _issue(
    severity: Severity,
    code: str,
    message: str,
    *,
    line: int | None = None,
    record_id: str | None = None,
) -> Issue:
    return Issue(severity, code, message, line, record_id)


def _inspect_file(path: Path) -> tuple[bool, str]:
    """Return gzip magic status and the SHA-256 of exact source bytes."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        first_bytes = handle.read(2)
        digest.update(first_bytes)
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return first_bytes == b"\x1f\x8b", digest.hexdigest()


def _iter_limited_lines(handle, max_line_length: int):
    """Yield text lines without allocating more than one permitted line."""

    while True:
        raw_line = handle.readline(max_line_length + 1)
        if not raw_line:
            return
        if len(raw_line) > max_line_length and not raw_line.endswith("\n"):
            raise ResourceLimitExceeded(
                "line_length_exceeded",
                f"Input line exceeds the {max_line_length}-character limit.",
            )
        yield raw_line


def _suffix_status(path: Path, is_gzip: bool) -> list[Issue]:
    suffixes = [suffix.lower() for suffix in path.suffixes]
    has_gzip_suffix = bool(suffixes and suffixes[-1] == ".gz")
    fasta_suffix = suffixes[-2] if has_gzip_suffix and len(suffixes) >= 2 else (suffixes[-1] if suffixes else "")
    issues: list[Issue] = []
    if has_gzip_suffix and not is_gzip:
        issues.append(
            _issue(
                "error",
                "gzip_suffix_not_gzip",
                "Filename ends in .gz but file bytes are not gzip.",
            )
        )
    if is_gzip and not has_gzip_suffix:
        issues.append(
            _issue(
                "warning",
                "gzip_without_suffix",
                "File bytes are gzip-compressed but filename does not end in .gz.",
            )
        )
    if fasta_suffix not in KNOWN_FASTA_SUFFIXES:
        issues.append(
            _issue(
                "warning",
                "unrecognised_fasta_suffix",
                "Filename does not use a recognised FASTA suffix; content will still be validated.",
            )
        )
    return issues


def _infer_alphabet(sequence: str) -> str:
    residues = frozenset(sequence.upper())
    if not residues:
        return "unknown"
    nucleotide_compatible = residues <= DNA_IUPAC or residues <= RNA_IUPAC
    protein_compatible = residues <= PROTEIN_IUPAC
    if nucleotide_compatible and protein_compatible:
        return "nucleotide_or_protein_ambiguous"
    if residues <= DNA_IUPAC:
        return "dna"
    if residues <= RNA_IUPAC:
        return "rna"
    if protein_compatible:
        return "protein"
    return "invalid"


def _allowed_residues(alphabet: Alphabet) -> frozenset[str]:
    return {
        "dna": DNA_IUPAC,
        "rna": RNA_IUPAC,
        "protein": PROTEIN_IUPAC,
    }[alphabet]  # type: ignore[index]


def _validate_sequence(
    sequence: str,
    *,
    line: int,
    record_id: str,
    options: LintOptions,
) -> tuple[str, list[Issue]]:
    normalized_for_validation = sequence.upper()
    inferred = _infer_alphabet(normalized_for_validation)
    issues: list[Issue] = []
    if options.alphabet != "auto":
        invalid = sorted(set(normalized_for_validation) - _allowed_residues(options.alphabet))
        if invalid:
            issues.append(
                _issue(
                    "error",
                    "alphabet_mismatch",
                    (
                        "Sequence contains residues incompatible with required "
                        f"{options.alphabet.upper()} alphabet: {''.join(invalid)}."
                    ),
                    line=line,
                    record_id=record_id,
                )
            )
    elif inferred == "invalid":
        invalid = sorted(set(normalized_for_validation) - PROTEIN_IUPAC)
        issues.append(
            _issue(
                "error",
                "invalid_residue",
                f"Sequence contains unsupported residues: {''.join(invalid)}.",
                line=line,
                record_id=record_id,
            )
        )

    ambiguity_symbols = (
        PROTEIN_AMBIGUITY
        if options.alphabet == "protein"
        else NUCLEOTIDE_AMBIGUITY
        if options.alphabet in {"dna", "rna"}
        else frozenset()
    )
    ambiguity = set(normalized_for_validation) & ambiguity_symbols
    if ambiguity:
        severity: Severity = "error" if options.disallow_ambiguity else "warning"
        issues.append(
            _issue(
                severity,
                "ambiguous_residue",
                f"Sequence contains ambiguous residues: {''.join(sorted(ambiguity))}.",
                line=line,
                record_id=record_id,
            )
        )
    return inferred, issues


def lint_fasta(path: str | Path, options: LintOptions | None = None) -> dict:
    """Validate a FASTA file and return a JSON-serialisable, stable report.

    This function never normalises residues or identifiers. I/O failures are
    returned as errors so callers can serialize one report for all outcomes.
    """

    options = options or LintOptions()
    source = Path(path)
    issues: list[Issue] = []
    records: list[RecordSummary] = []
    input_sha256: str | None = None
    detected_header_count = 0
    try:
        if not source.is_file():
            raise FileNotFoundError(source)
        if source.stat().st_size > options.max_input_bytes:
            raise ResourceLimitExceeded(
                "input_size_exceeded",
                f"Input file exceeds the {options.max_input_bytes}-byte limit.",
            )
        is_gzip, input_sha256 = _inspect_file(source)
        issues.extend(_suffix_status(source, is_gzip))
        opener = gzip.open if is_gzip else open
        with opener(source, "rt", encoding="utf-8", errors="strict", newline=None) as handle:
            current_id: str | None = None
            current_description = ""
            current_has_content_line = False
            current_sequence_length = 0
            current_sequence_digest = hashlib.sha256()
            header_line: int | None = None
            seen_ids: set[str] = set()
            preflight_records: list[PreflightRecord] = []
            total_sequence_length = 0

            def finish_record() -> None:
                nonlocal current_id, current_description, current_has_content_line
                nonlocal current_sequence_length, current_sequence_digest, header_line
                if current_id is None:
                    return
                if not current_has_content_line:
                    issues.append(
                        _issue(
                            "error",
                            "empty_sequence",
                            "FASTA record has no sequence residues.",
                            line=header_line,
                            record_id=current_id,
                        )
                    )
                    return
                if not current_sequence_length:
                    return
                preflight_records.append(
                    PreflightRecord(
                        record_id=current_id,
                        description=current_description,
                        header_line=header_line or 0,
                        sequence_length=current_sequence_length,
                        sequence_sha256=current_sequence_digest.hexdigest(),
                    )
                )

            for line_number, raw_line in enumerate(
                _iter_limited_lines(handle, options.max_line_length),
                start=1,
            ):
                line = raw_line.rstrip("\r\n")
                if not line:
                    if current_id is None:
                        issues.append(
                            _issue(
                                "error",
                                "blank_line_before_header",
                                "Blank line appears before the first FASTA header.",
                                line=line_number,
                            )
                        )
                        continue
                    issues.append(
                        _issue(
                            "warning",
                            "blank_line",
                            "Blank line found in FASTA input.",
                            line=line_number,
                            record_id=current_id,
                        )
                    )
                    continue
                if line.startswith(">"):
                    finish_record()
                    detected_header_count += 1
                    if detected_header_count > options.max_records:
                        raise ResourceLimitExceeded(
                            "record_limit_exceeded",
                            f"Input exceeds the {options.max_records}-record limit.",
                        )
                    header = line[1:]
                    if not header.strip():
                        issues.append(
                            _issue(
                                "error",
                                "empty_header",
                                "FASTA header is empty.",
                                line=line_number,
                            )
                        )
                        current_id = f"__invalid_header_line_{line_number}"
                        current_description = ""
                    else:
                        parts = header.split(maxsplit=1)
                        current_id = parts[0]
                        current_description = parts[1] if len(parts) == 2 else ""
                        if current_id in seen_ids:
                            issues.append(
                                _issue(
                                    "error",
                                    "duplicate_identifier",
                                    "FASTA identifier occurs more than once.",
                                    line=line_number,
                                    record_id=current_id,
                                )
                            )
                        seen_ids.add(current_id)
                    current_has_content_line = False
                    current_sequence_length = 0
                    current_sequence_digest = hashlib.sha256()
                    header_line = line_number
                    continue
                if current_id is None:
                    issues.append(
                        _issue(
                            "error",
                            "sequence_before_header",
                            "Sequence data appears before the first FASTA header.",
                            line=line_number,
                        )
                    )
                    continue
                if any(character.isspace() for character in line):
                    current_has_content_line = True
                    issues.append(
                        _issue(
                            "error",
                            "whitespace_in_sequence",
                            "Sequence line contains whitespace and was not normalised.",
                            line=line_number,
                            record_id=current_id,
                        )
                    )
                else:
                    current_has_content_line = True
                    current_sequence_length += len(line)
                    total_sequence_length += len(line)
                    if current_sequence_length > options.max_sequence_length:
                        raise ResourceLimitExceeded(
                            "sequence_length_exceeded",
                            f"Record exceeds the {options.max_sequence_length}-residue limit.",
                        )
                    if total_sequence_length > options.max_total_sequence_length:
                        raise ResourceLimitExceeded(
                            "total_sequence_length_exceeded",
                            f"Input exceeds the {options.max_total_sequence_length}-residue limit.",
                        )
                    current_sequence_digest.update(line.encode("utf-8"))
            finish_record()

        if detected_header_count == 0:
            issues.append(_issue("error", "no_records", "FASTA input contains no records."))

        # Biopython is the authoritative FASTA record parser. The preceding
        # byte-preserving pass supplies policy diagnostics (for example,
        # whitespace inside a sequence) that permissive FASTA readers may
        # otherwise normalise silently.
        suppress_record_codes = {
            "duplicate_identifier",
            "empty_header",
            "empty_sequence",
            "sequence_before_header",
            "whitespace_in_sequence",
            "blank_line_before_header",
            "gzip_suffix_not_gzip",
        }
        has_suppressing_errors = any(
            issue.severity == "error" and issue.code in suppress_record_codes for issue in issues
        )
        if has_suppressing_errors:
            # A partial summary can be mistaken for a valid downstream input.
            # Structural errors therefore suppress all record-level output.
            preflight_records = []
        else:
            with opener(source, "rt", encoding="utf-8", errors="strict", newline=None) as handle:
                parser = iter(SeqIO.parse(handle, "fasta"))
                parser_integrity_failed = False
                for expected in preflight_records:
                    record = next(parser, None)
                    if record is None:
                        issues.append(
                            _issue(
                                "error",
                                "parser_record_count_mismatch",
                                "Biopython ended before the preflight record count.",
                            )
                        )
                        parser_integrity_failed = True
                        break
                    sequence = str(record.seq)
                    sequence_sha256 = hashlib.sha256(sequence.encode("utf-8")).hexdigest()
                    if (
                        record.id != expected.record_id
                        or len(sequence) != expected.sequence_length
                        or sequence_sha256 != expected.sequence_sha256
                    ):
                        issues.append(
                            _issue(
                                "error",
                                "parser_record_mismatch",
                                "Biopython record content differs from the preflight record.",
                                line=expected.header_line,
                                record_id=expected.record_id,
                            )
                        )
                        parser_integrity_failed = True
                        break
                    inferred, sequence_issues = _validate_sequence(
                        sequence,
                        line=expected.header_line,
                        record_id=expected.record_id,
                        options=options,
                    )
                    issues.extend(sequence_issues)
                    records.append(
                        RecordSummary(
                            record_id=expected.record_id,
                            description=expected.description,
                            sequence_length=len(sequence),
                            inferred_alphabet=inferred,
                            sha256=sequence_sha256,
                        )
                    )
                if not parser_integrity_failed and next(parser, None) is not None:
                    issues.append(
                        _issue(
                            "error",
                            "parser_record_count_mismatch",
                            "Biopython produced more records than the preflight record count.",
                        )
                    )

        integrity_issue_codes = {"parser_record_count_mismatch", "parser_record_mismatch"}
        if any(issue.code in integrity_issue_codes for issue in issues):
            # A parser disagreement invalidates the entire input boundary; a
            # prefix of apparently valid records must never escape downstream.
            records.clear()
    except ResourceLimitExceeded as error:
        issues.append(_issue("error", error.code, str(error)))
        records.clear()
    except (OSError, UnicodeError, EOFError, ValueError) as error:
        issues.append(_issue("error", "unreadable_input", f"Could not read input: {error}."))

    error_count = sum(issue.severity == "error" for issue in issues)
    warning_count = sum(issue.severity == "warning" for issue in issues)
    return {
        "schema_version": "1.0",
        "tool": "humaninfinder",
        "report_type": "fasta_lint",
        "created_at_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "input": {
            "path": str(source),
            "sha256": input_sha256,
            "requested_alphabet": options.alphabet,
            "disallow_ambiguity": options.disallow_ambiguity,
            "limits": {
                "max_input_bytes": options.max_input_bytes,
                "max_records": options.max_records,
                "max_line_length": options.max_line_length,
                "max_sequence_length": options.max_sequence_length,
                "max_total_sequence_length": options.max_total_sequence_length,
            },
        },
        "summary": {
            "record_count": len(records),
            "detected_header_count": detected_header_count,
            "error_count": error_count,
            "warning_count": warning_count,
            "valid": error_count == 0,
        },
        "records": [asdict(record) for record in records],
        "issues": [asdict(issue) for issue in issues],
    }
