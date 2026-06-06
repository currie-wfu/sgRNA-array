"""5' array barcode generation and validation.

Each array gets a single 12-nt barcode upstream of the position-1 stem-I. The barcode
is for bulk PCR/Sanger verification of array identity from gDNA; it is NOT a 10x
direct-capture handle (which lives in the scaffold's Hairpin-1 loop when CS1 is enabled).

Generated barcodes satisfy:
- Exact length
- ACGT alphabet
- GC fraction in [DEFAULT_BARCODE_GC_MIN, DEFAULT_BARCODE_GC_MAX]
- No homopolymer run > DEFAULT_BARCODE_MAX_HOMOPOLYMER
- No PaqCI recognition site
- No CS1 / CS2 homology (Hamming distance threshold within sliding window)
- Minimum pairwise Hamming distance with other barcodes in the same batch
"""

from __future__ import annotations

import random

from sgrna_array.constants import (
    CS1_INSGRNA,
    CS2_INSGRNA,
    DEFAULT_BARCODE_GC_MAX,
    DEFAULT_BARCODE_GC_MIN,
    DEFAULT_BARCODE_LEN,
    DEFAULT_BARCODE_MAX_HOMOPOLYMER,
    DEFAULT_BARCODE_MIN_HAMMING,
)
from sgrna_array.enzymes import scan_paqci_sites
from sgrna_array.validator import Severity, ValidationResult


def hamming(a: str, b: str) -> int:
    """Hamming distance between two equal-length strings (case-sensitive)."""
    if len(a) != len(b):
        raise ValueError(f"Hamming distance requires equal-length strings; got {len(a)} vs {len(b)}")
    return sum(1 for x, y in zip(a, b, strict=True) if x != y)


def generate_barcode(
    rng: random.Random,
    length: int = DEFAULT_BARCODE_LEN,
    existing: list[str] | None = None,
    max_attempts: int = 10_000,
) -> str:
    """Generate a single barcode satisfying all constraints.

    Args:
        rng: A `random.Random` instance. Pass a seeded one for reproducibility.
        length: Barcode length in nt.
        existing: Other barcodes already issued in this batch; new barcode must be
            ≥ DEFAULT_BARCODE_MIN_HAMMING from each.
        max_attempts: Maximum number of random draws before giving up.

    Returns:
        A new barcode string (uppercase ACGT).

    Raises:
        RuntimeError: If no valid barcode is found in `max_attempts` tries.
    """
    existing = existing or []
    bases = "ACGT"
    for _ in range(max_attempts):
        candidate = "".join(rng.choice(bases) for _ in range(length))
        if _passes_intrinsic(candidate) and _passes_collisions(candidate, existing):
            return candidate
    raise RuntimeError(
        f"Failed to generate a valid barcode in {max_attempts} attempts; "
        "loosen constraints or check for over-constrained inputs"
    )


def validate_barcode(seq: str, existing: list[str] | None = None) -> list[ValidationResult]:
    """Validate a user-supplied barcode and return any findings.

    Same constraints as `generate_barcode` but reported as validation results instead of
    enforced by rejection sampling.
    """
    existing = existing or []
    results: list[ValidationResult] = []
    if len(seq) != DEFAULT_BARCODE_LEN:
        results.append(
            ValidationResult(
                severity=Severity.WARNING,
                code="barcode_length",
                message=(
                    f"Barcode length {len(seq)} differs from default {DEFAULT_BARCODE_LEN}"
                ),
            )
        )
    gc = _gc_fraction(seq)
    if not (DEFAULT_BARCODE_GC_MIN <= gc <= DEFAULT_BARCODE_GC_MAX):
        results.append(
            ValidationResult(
                severity=Severity.WARNING,
                code="barcode_gc",
                message=(
                    f"Barcode GC {gc:.0%} outside recommended "
                    f"[{DEFAULT_BARCODE_GC_MIN:.0%}, {DEFAULT_BARCODE_GC_MAX:.0%}]"
                ),
            )
        )
    run = _longest_homopolymer(seq)
    if run > DEFAULT_BARCODE_MAX_HOMOPOLYMER:
        results.append(
            ValidationResult(
                severity=Severity.WARNING,
                code="barcode_homopolymer",
                message=(
                    f"Barcode contains {run}-nt homopolymer run "
                    f"(threshold {DEFAULT_BARCODE_MAX_HOMOPOLYMER})"
                ),
            )
        )
    if scan_paqci_sites(seq):
        results.append(
            ValidationResult(
                severity=Severity.ERROR,
                code="barcode_paqci_site",
                message="Barcode contains internal PaqCI recognition site",
            )
        )
    if _has_cs_homology(seq):
        results.append(
            ValidationResult(
                severity=Severity.WARNING,
                code="barcode_cs_homology",
                message="Barcode has high homology to CS1 or CS2 capture sequence",
            )
        )
    for other in existing:
        if len(other) == len(seq) and hamming(seq, other) < DEFAULT_BARCODE_MIN_HAMMING:
            results.append(
                ValidationResult(
                    severity=Severity.WARNING,
                    code="barcode_collision",
                    message=(
                        f"Barcode is within Hamming distance "
                        f"{hamming(seq, other)} of existing barcode {other!r}"
                    ),
                )
            )
            break
    return results


def _passes_intrinsic(seq: str) -> bool:
    """Run the per-barcode constraints (GC, homopolymer, no PaqCI, no CS homology)."""
    if not (DEFAULT_BARCODE_GC_MIN <= _gc_fraction(seq) <= DEFAULT_BARCODE_GC_MAX):
        return False
    if _longest_homopolymer(seq) > DEFAULT_BARCODE_MAX_HOMOPOLYMER:
        return False
    if scan_paqci_sites(seq):
        return False
    if _has_cs_homology(seq):
        return False
    return True


def _passes_collisions(seq: str, existing: list[str]) -> bool:
    return all(
        hamming(seq, other) >= DEFAULT_BARCODE_MIN_HAMMING
        for other in existing
        if len(other) == len(seq)
    )


def _gc_fraction(seq: str) -> float:
    return sum(1 for c in seq.upper() if c in "GC") / len(seq)


def _longest_homopolymer(seq: str) -> int:
    longest = 1
    current = 1
    for i in range(1, len(seq)):
        if seq[i].upper() == seq[i - 1].upper():
            current += 1
            longest = max(longest, current)
        else:
            current = 1
    return longest


def _has_cs_homology(seq: str, threshold: int = 4) -> bool:
    """Check if `seq` has a sliding window within Hamming distance `threshold` of CS1/CS2.

    Compares every length-(len(seq)) window of CS1 and CS2 against `seq`.
    """
    seq_upper = seq.upper()
    for cs in (CS1_INSGRNA, CS2_INSGRNA):
        cs_upper = cs.upper()
        if len(seq_upper) > len(cs_upper):
            continue
        for i in range(len(cs_upper) - len(seq_upper) + 1):
            window = cs_upper[i : i + len(seq_upper)]
            if hamming(seq_upper, window) < threshold:
                return True
    return False
