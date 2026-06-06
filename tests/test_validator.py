"""Validator tests: length / alphabet / PaqCI / GC / homopolymer rules."""

from __future__ import annotations

from sgrna_array.validator import Severity, has_errors, validate_crrna


def test_valid_crrna_produces_no_findings() -> None:
    results = validate_crrna("ACGTACGTACGTACGTACGT")
    assert results == []


def test_wrong_length_is_hard_error() -> None:
    results = validate_crrna("ACGT")
    assert has_errors(results)
    assert any(r.code == "bad_length" for r in results)


def test_invalid_alphabet_is_hard_error() -> None:
    results = validate_crrna("ACGTACGTACGTACGTACGN")  # contains N
    assert has_errors(results)
    assert any(r.code == "bad_alphabet" for r in results)


def test_internal_paqci_fwd_is_hard_error() -> None:
    # `CACCTGC` embedded in an otherwise valid 20-nt spacer.
    seq = "AAACACCTGCAAATTACGTT"  # 20 nt, contains CACCTGC
    assert len(seq) == 20
    results = validate_crrna(seq)
    assert has_errors(results)
    assert any(r.code == "internal_paqci_site" for r in results)


def test_internal_paqci_rev_is_hard_error() -> None:
    seq = "AAAGCAGGTGAAATTACGTT"  # 20 nt, contains GCAGGTG
    assert len(seq) == 20
    results = validate_crrna(seq)
    assert has_errors(results)
    assert any(r.code == "internal_paqci_site" for r in results)


def test_low_gc_warns_but_does_not_fail() -> None:
    results = validate_crrna("ATATATATATATATATATAT")  # 0% GC
    assert not has_errors(results)
    assert any(r.code == "low_gc" and r.severity is Severity.WARNING for r in results)


def test_high_gc_warns_but_does_not_fail() -> None:
    results = validate_crrna("GCGCGCGCGCGCGCGCGCGC")  # 100% GC
    assert not has_errors(results)
    assert any(r.code == "high_gc" and r.severity is Severity.WARNING for r in results)


def test_long_homopolymer_warns_but_does_not_fail() -> None:
    results = validate_crrna("AAAAAAACGTACGTACGTAC")  # 7-run of A
    assert not has_errors(results)
    assert any(r.code == "long_homopolymer" and r.severity is Severity.WARNING for r in results)


def test_validator_skips_downstream_checks_when_basics_fail() -> None:
    """If the crRNA fails length or alphabet, downstream checks should not run."""
    results = validate_crrna("ACGT")  # bad length
    # Should only have the bad_length error, not GC or homopolymer warnings.
    codes = {r.code for r in results}
    assert codes == {"bad_length"}
