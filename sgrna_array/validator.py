"""crRNA validation: length, alphabet, internal PaqCI sites, GC content, homopolymers.

Splice-site scanning is deferred to v1.1 (see roadmap). The v1 validator catches the
hard-fail issues that would break assembly or synthesis.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from sgrna_array.constants import (
    CRRNA_GC_WARN_MAX,
    CRRNA_GC_WARN_MIN,
    CRRNA_HOMOPOLYMER_WARN_THRESHOLD,
    CRRNA_LEN,
    DNA_ALPHABET,
)
from sgrna_array.enzymes import scan_paqci_sites


class Severity(Enum):
    """Validation result severity. Hard fails block assembly; warnings are reported."""

    ERROR = "error"
    WARNING = "warning"


@dataclass(frozen=True)
class ValidationResult:
    """A single validation finding against an input crRNA or full fragment."""

    severity: Severity
    code: str
    message: str

    def __str__(self) -> str:  # pragma: no cover - trivial
        return f"[{self.severity.value.upper()}] {self.code}: {self.message}"


def validate_crrna(seq: str) -> list[ValidationResult]:
    """Run all v1 crRNA validations and return the (possibly empty) list of findings."""
    results: list[ValidationResult] = []
    results.extend(_check_length(seq))
    results.extend(_check_alphabet(seq))
    # Only run downstream checks if the basics passed, since they assume well-formed input.
    if not any(r.severity is Severity.ERROR for r in results):
        results.extend(_check_no_internal_paqci(seq))
        results.extend(_check_gc_content(seq))
        results.extend(_check_homopolymer(seq))
    return results


def has_errors(results: list[ValidationResult]) -> bool:
    """Return True if any finding is a hard error."""
    return any(r.severity is Severity.ERROR for r in results)


def _check_length(seq: str) -> list[ValidationResult]:
    if len(seq) != CRRNA_LEN:
        return [
            ValidationResult(
                severity=Severity.ERROR,
                code="bad_length",
                message=f"crRNA must be exactly {CRRNA_LEN} nt; got {len(seq)}",
            )
        ]
    return []


def _check_alphabet(seq: str) -> list[ValidationResult]:
    invalid = sorted({ch for ch in seq if ch not in DNA_ALPHABET})
    if invalid:
        return [
            ValidationResult(
                severity=Severity.ERROR,
                code="bad_alphabet",
                message=f"crRNA contains non-ACGT characters: {invalid!r}",
            )
        ]
    return []


def _check_no_internal_paqci(seq: str) -> list[ValidationResult]:
    hits = scan_paqci_sites(seq)
    if hits:
        return [
            ValidationResult(
                severity=Severity.ERROR,
                code="internal_paqci_site",
                message=(
                    "crRNA contains internal PaqCI recognition site(s) that would break "
                    f"Golden Gate assembly: {hits}"
                ),
            )
        ]
    return []


def _check_gc_content(seq: str) -> list[ValidationResult]:
    gc = sum(1 for c in seq.upper() if c in "GC") / len(seq)
    if gc < CRRNA_GC_WARN_MIN:
        return [
            ValidationResult(
                severity=Severity.WARNING,
                code="low_gc",
                message=f"crRNA GC content {gc:.0%} below recommended {CRRNA_GC_WARN_MIN:.0%}",
            )
        ]
    if gc > CRRNA_GC_WARN_MAX:
        return [
            ValidationResult(
                severity=Severity.WARNING,
                code="high_gc",
                message=f"crRNA GC content {gc:.0%} above recommended {CRRNA_GC_WARN_MAX:.0%}",
            )
        ]
    return []


def _check_homopolymer(seq: str) -> list[ValidationResult]:
    """Soft warn on homopolymer runs ≥ CRRNA_HOMOPOLYMER_WARN_THRESHOLD."""
    longest_run = 1
    longest_base = seq[0]
    current_run = 1
    for i in range(1, len(seq)):
        if seq[i].upper() == seq[i - 1].upper():
            current_run += 1
            if current_run > longest_run:
                longest_run = current_run
                longest_base = seq[i]
        else:
            current_run = 1
    if longest_run >= CRRNA_HOMOPOLYMER_WARN_THRESHOLD:
        return [
            ValidationResult(
                severity=Severity.WARNING,
                code="long_homopolymer",
                message=(
                    f"crRNA contains a {longest_run}-nt run of {longest_base.upper()} "
                    f"(threshold ≥ {CRRNA_HOMOPOLYMER_WARN_THRESHOLD}); "
                    "may affect synthesis quality"
                ),
            )
        ]
    return []
