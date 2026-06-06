"""Barcode generation and validation tests."""

from __future__ import annotations

import random

import pytest

from sgrna_array.barcode import generate_barcode, hamming, validate_barcode
from sgrna_array.constants import DEFAULT_BARCODE_LEN, DEFAULT_BARCODE_MIN_HAMMING


def test_hamming_basics() -> None:
    assert hamming("ACGT", "ACGT") == 0
    assert hamming("ACGT", "ACGA") == 1
    assert hamming("ACGT", "TGCA") == 4


def test_hamming_requires_equal_length() -> None:
    with pytest.raises(ValueError, match="equal-length"):
        hamming("ACGT", "ACG")


def test_generate_barcode_is_deterministic_with_seed() -> None:
    rng1 = random.Random(42)
    rng2 = random.Random(42)
    assert generate_barcode(rng1) == generate_barcode(rng2)


def test_generate_barcode_default_length() -> None:
    bc = generate_barcode(random.Random(0))
    assert len(bc) == DEFAULT_BARCODE_LEN


def test_generate_barcode_satisfies_constraints() -> None:
    """Generated barcode must pass `validate_barcode` with no warnings except possibly
    collision warnings (but here we pass no existing barcodes)."""
    bc = generate_barcode(random.Random(123))
    findings = validate_barcode(bc)
    # No warnings or errors expected for a freshly generated barcode against empty batch.
    assert findings == []


def test_generate_barcode_respects_existing_collisions() -> None:
    """A new barcode must be at least `DEFAULT_BARCODE_MIN_HAMMING` from each existing one."""
    rng = random.Random(7)
    batch: list[str] = []
    for _ in range(5):
        bc = generate_barcode(rng, existing=batch)
        for existing in batch:
            assert hamming(bc, existing) >= DEFAULT_BARCODE_MIN_HAMMING
        batch.append(bc)


def test_validate_barcode_flags_paqci_site() -> None:
    bc = "CACCTGCACGTA"  # 12 nt, embeds CACCTGC
    findings = validate_barcode(bc)
    assert any(f.code == "barcode_paqci_site" for f in findings)


def test_validate_barcode_flags_long_homopolymer() -> None:
    bc = "AAAAACGTACGT"  # 5-run of A
    findings = validate_barcode(bc)
    assert any(f.code == "barcode_homopolymer" for f in findings)


def test_validate_barcode_flags_collision() -> None:
    findings = validate_barcode("ACGTACGTACGT", existing=["ACGTACGTACGA"])  # Hamming = 1
    assert any(f.code == "barcode_collision" for f in findings)
