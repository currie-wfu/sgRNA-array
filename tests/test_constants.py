"""Smoke tests that every locked constant has the documented length / shape.

These tests are the first line of defense against accidental sequence corruption when
editing constants.py. They don't validate scientific correctness (that's the dossier's
job) — they just catch obvious typos like a deleted nucleotide.
"""

from __future__ import annotations

import re

from sgrna_array import constants


def test_paqci_recognition_sites_are_reverse_complements() -> None:
    table = str.maketrans("ACGT", "TGCA")
    assert (
        constants.PAQCI_RECOGNITION_FWD.translate(table)[::-1]
        == constants.PAQCI_RECOGNITION_REV
    )


def test_hh_catalytic_core_length() -> None:
    assert len(constants.HH_CATALYTIC_CORE) == 37


def test_hdv_ribozyme_length() -> None:
    assert len(constants.HDV_RIBOZYME) == 68


def test_dna_only_alphabet_for_sequence_constants() -> None:
    """Every sequence constant should consist of A/C/G/T only (case-insensitive)."""
    seq_attrs = [
        "HH_CATALYTIC_CORE",
        "HDV_RIBOZYME",
        "SGRNA_SCAFFOLD_DEWEIRDT",
        "SGRNA_SCAFFOLD_DEWEIRDT_CS1",
        "SGRNA_SCAFFOLD_HSU_MALI",
        "CS1_INSGRNA",
        "CS1_GELBEAD",
        "CS2_INSGRNA",
        "CS2_GELBEAD",
        "CS_HAIRPIN_STEM",
        "B5_OVERHANG",
        "B3_OVERHANG",
        "STITCHER_SPACER",
        "PAQCI_RECOGNITION_FWD",
        "PAQCI_RECOGNITION_REV",
        "PAQCI_RECOGNITION_TO_OVERHANG_SPACER",
    ]
    for attr in seq_attrs:
        value = getattr(constants, attr)
        assert re.fullmatch(r"[ACGTacgt]+", value), f"{attr} contains non-ACGT chars: {value!r}"


def test_deweirdt_scaffold_length_and_mutations() -> None:
    s = constants.SGRNA_SCAFFOLD_DEWEIRDT
    assert len(s) == 76
    # Position 5 should be G (DeWeirdt T5→G mutation; 0-indexed position 4).
    assert s[4].lower() == "g", f"DeWeirdt scaffold position 5 should be G, got {s[4]!r}"
    # Position 26 should be C (DeWeirdt A26→C compensatory mutation; 0-indexed position 25).
    assert s[25].lower() == "c", f"DeWeirdt scaffold position 26 should be C, got {s[25]!r}"


def test_hsu_mali_scaffold_length_and_no_deweirdt_mutations() -> None:
    s = constants.SGRNA_SCAFFOLD_HSU_MALI
    assert len(s) == 76
    # Position 5 should be T (the Pol III termination signal residue).
    assert s[4].lower() == "t"
    # Position 26 should be A (the pre-DeWeirdt baseline).
    assert s[25].lower() == "a"


def test_cs1_modified_scaffold_contains_cs1_in_hairpin() -> None:
    s = constants.SGRNA_SCAFFOLD_DEWEIRDT_CS1
    assert len(s) == 102, f"CS1-modified scaffold should be 102 nt, got {len(s)}"
    assert constants.CS1_INSGRNA.lower() in s.lower()
    # Flanked by the ggcc-...-ggcc stem.
    cs1_stem = constants.CS_HAIRPIN_STEM + constants.CS1_INSGRNA + constants.CS_HAIRPIN_STEM
    assert cs1_stem.lower() in s.lower()


def test_capture_sequences_are_reverse_complements() -> None:
    table = str.maketrans("ACGT", "TGCA")
    assert constants.CS1_INSGRNA.translate(table)[::-1] == constants.CS1_GELBEAD
    assert constants.CS2_INSGRNA.translate(table)[::-1] == constants.CS2_GELBEAD


def test_overhang_set_size_and_uniqueness() -> None:
    all_overhangs = [constants.B5_OVERHANG, *constants.J_OVERHANGS, constants.B3_OVERHANG]
    assert len(all_overhangs) == 13, "Expect 13 junctions (B5 + 11 J + B3)"
    assert len(set(all_overhangs)) == 13, "Overhang set must be unique"
    for oh in all_overhangs:
        assert len(oh) == 4 and re.fullmatch(r"[ACGT]+", oh), f"Bad overhang: {oh!r}"


def test_b5_b3_are_vector_locked_values() -> None:
    assert constants.B5_OVERHANG == "ACGG"
    assert constants.B3_OVERHANG == "GAGC"


def test_stitcher_junction_map_uses_valid_j_overhangs() -> None:
    j_set = set(constants.J_OVERHANGS)
    for size, overhang in constants.STITCHER_JUNCTION_BY_STOP_SIZE.items():
        assert size in {4, 6, 8, 10}
        assert overhang in j_set, f"Stitcher overhang {overhang!r} for size {size} not in J set"


def test_supported_array_sizes() -> None:
    assert constants.SUPPORTED_ARRAY_SIZES == (1, 2, 4, 6, 8, 10, 12)


def test_terminal_uses_b3_directly_covers_endpoint_and_singletons() -> None:
    # 1 and 2: too short to use library stitchers; terminal fragment goes direct to B3.
    # 12: full-library natural endpoint at B3.
    assert constants.TERMINAL_USES_B3_DIRECTLY == frozenset({1, 2, 12})
    # No overlap with the stitcher map (each size belongs to exactly one regime).
    assert not (
        constants.TERMINAL_USES_B3_DIRECTLY
        & constants.STITCHER_JUNCTION_BY_STOP_SIZE.keys()
    )
    # All supported sizes are accounted for by one regime or the other.
    assert (
        constants.TERMINAL_USES_B3_DIRECTLY
        | constants.STITCHER_JUNCTION_BY_STOP_SIZE.keys()
        == set(constants.SUPPORTED_ARRAY_SIZES)
    )
