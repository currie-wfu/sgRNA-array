"""Golden-fragment tests for the assembler.

Verifies the end-to-end `build_fragment` and `build_array` flow against deterministic
expected outputs computed from the locked constants. Uses the Shh gRNA#3 example crRNA
from the host vector (`ggtgacgcgtgtgtacctgg`) since its stem-I (`gtcacc`) is dossier-verified.
"""

from __future__ import annotations

import pytest

from sgrna_array import constants
from sgrna_array.assembler import (
    build_array,
    build_fragment,
    build_stitcher,
    junction_overhang,
    left_overhang,
)
from sgrna_array.ribozyme import build_hh_stem, reverse_complement

# Reference inputs.
EXAMPLE_CRRNA = "ggtgacgcgtgtgtacctgg"
EXPECTED_STEM_I = "gtcacc"


# ---------------------------------------------------------------------------
# Ribozyme helpers
# ---------------------------------------------------------------------------


def test_reverse_complement_roundtrip() -> None:
    assert reverse_complement("ACGT") == "ACGT"  # palindrome
    assert reverse_complement(EXAMPLE_CRRNA) == "ccaggtacacacgcgtcacc"
    # Double-RC = identity
    assert reverse_complement(reverse_complement(EXAMPLE_CRRNA)) == EXAMPLE_CRRNA


def test_hh_stem_matches_dossier_for_shh_grna() -> None:
    assert build_hh_stem(EXAMPLE_CRRNA) == EXPECTED_STEM_I


def test_hh_stem_rejects_too_short_input() -> None:
    with pytest.raises(ValueError, match="at least"):
        build_hh_stem("ACGTA")  # 5 nt < HH_STEM_I_LEN


# ---------------------------------------------------------------------------
# Overhang junction logic
# ---------------------------------------------------------------------------


def test_left_overhang_position_1_is_b5() -> None:
    assert left_overhang(1) == constants.B5_OVERHANG


def test_left_overhang_inner_positions_use_j_overhangs() -> None:
    # Position N's left overhang = J(N-1), which is index N-2 in the tuple.
    assert left_overhang(2) == constants.J_OVERHANGS[0]
    assert left_overhang(7) == constants.J_OVERHANGS[5]
    assert left_overhang(12) == constants.J_OVERHANGS[10]


def test_junction_overhang_terminal_position_12_is_b3() -> None:
    assert junction_overhang(12, 12) == constants.B3_OVERHANG


def test_junction_overhang_terminal_position_for_sizes_1_and_2_is_b3() -> None:
    # Sizes 1 and 2: terminal fragment ligates straight to backbone B3.
    assert junction_overhang(1, 1) == constants.B3_OVERHANG
    assert junction_overhang(2, 2) == constants.B3_OVERHANG


def test_junction_overhang_terminal_position_for_sub12_arrays() -> None:
    # Terminal of a stop-at-4 array goes to the J4 overhang (= stitcher 5' overhang).
    assert junction_overhang(4, 4) == constants.STITCHER_JUNCTION_BY_STOP_SIZE[4]
    assert junction_overhang(6, 6) == constants.STITCHER_JUNCTION_BY_STOP_SIZE[6]


def test_junction_overhang_internal_position() -> None:
    # Inside an array, the 3' overhang of position N = J_N (index N-1).
    assert junction_overhang(1, 12) == constants.J_OVERHANGS[0]
    assert junction_overhang(5, 12) == constants.J_OVERHANGS[4]


# ---------------------------------------------------------------------------
# Single-fragment assembly
# ---------------------------------------------------------------------------


def test_build_fragment_pos1_array4_no_barcode_no_cs1() -> None:
    frag = build_fragment(
        crrna=EXAMPLE_CRRNA,
        position=1,
        array_size=4,
        gene_label="Shh",
    )
    assert frag.position == 1
    assert frag.gene_label == "Shh"
    assert frag.left_overhang == constants.B5_OVERHANG
    # Position 1 of an array_size-4 array → 3' overhang is J1.
    assert frag.right_overhang == constants.J_OVERHANGS[0]
    assert frag.barcode is None
    assert frag.use_cs1 is False
    # The RNA-coding core must contain stem-I, crRNA, scaffold, and HDV.
    assert EXPECTED_STEM_I in frag.ordered_dna.lower()
    assert EXAMPLE_CRRNA.lower() in frag.ordered_dna.lower()
    assert constants.SGRNA_SCAFFOLD_DEWEIRDT.lower() in frag.ordered_dna.lower()
    assert constants.HDV_RIBOZYME.lower() in frag.ordered_dna.lower()
    # The cassette must have the PaqCI recognition sites in canonical positions.
    assert constants.PAQCI_RECOGNITION_FWD in frag.ordered_dna
    assert constants.PAQCI_RECOGNITION_REV in frag.ordered_dna


def test_build_fragment_length_is_predictable() -> None:
    frag = build_fragment(crrna=EXAMPLE_CRRNA, position=1, array_size=12)
    # Expected length: 5 pad + 7 PaqCI + 1 spacer + 4 overhang + 6 stem-I + 37 HH core
    # + 20 crRNA + 76 scaffold + 68 HDV + 4 overhang + 1 spacer + 7 PaqCI + 5 pad
    # = 241 nt
    expected_len = (
        constants.PAQCI_FLANKING_PAD_LEN
        + len(constants.PAQCI_RECOGNITION_FWD)
        + len(constants.PAQCI_RECOGNITION_TO_OVERHANG_SPACER)
        + 4  # left overhang
        + constants.HH_STEM_I_LEN
        + len(constants.HH_CATALYTIC_CORE)
        + constants.CRRNA_LEN
        + len(constants.SGRNA_SCAFFOLD_DEWEIRDT)
        + len(constants.HDV_RIBOZYME)
        + 4  # right overhang
        + len(constants.PAQCI_RECOGNITION_TO_OVERHANG_SPACER)
        + len(constants.PAQCI_RECOGNITION_REV)
        + constants.PAQCI_FLANKING_PAD_LEN
    )
    assert len(frag.ordered_dna) == expected_len == 241


def test_build_fragment_with_barcode_includes_it_only_at_position_1() -> None:
    bc = "ACGTAGCATCGT"  # 12 nt, balanced GC, no homopolymer (just for the test)
    frag1 = build_fragment(
        crrna=EXAMPLE_CRRNA, position=1, array_size=4, barcode=bc
    )
    frag2 = build_fragment(
        crrna=EXAMPLE_CRRNA, position=2, array_size=4, barcode=bc
    )
    assert frag1.barcode == bc
    assert frag2.barcode is None  # silently dropped for non-position-1
    assert bc in frag1.ordered_dna
    assert bc not in frag2.ordered_dna
    assert len(frag1.ordered_dna) == len(frag2.ordered_dna) + len(bc)


def test_build_fragment_with_cs1_uses_modified_scaffold() -> None:
    frag = build_fragment(crrna=EXAMPLE_CRRNA, position=1, array_size=4, use_cs1=True)
    assert frag.use_cs1 is True
    assert constants.SGRNA_SCAFFOLD_DEWEIRDT_CS1.lower() in frag.ordered_dna.lower()
    assert constants.CS1_INSGRNA.lower() in frag.ordered_dna.lower()


def test_build_fragment_rejects_bad_position() -> None:
    with pytest.raises(ValueError, match="position"):
        build_fragment(crrna=EXAMPLE_CRRNA, position=0, array_size=4)
    with pytest.raises(ValueError, match="position"):
        build_fragment(crrna=EXAMPLE_CRRNA, position=5, array_size=4)


def test_build_fragment_rejects_unsupported_array_size() -> None:
    with pytest.raises(ValueError, match="array_size"):
        build_fragment(crrna=EXAMPLE_CRRNA, position=1, array_size=7)


# ---------------------------------------------------------------------------
# Array-level assembly + stitchers
# ---------------------------------------------------------------------------


def _dummy_crrnas(n: int) -> list[str]:
    """Generate n distinct 20-nt crRNAs for size tests. Sequences are biologically
    meaningless but pass the validator (no PaqCI sites, normal GC)."""
    return [
        "ACGTACGTAC" + base * 10
        for base in ("A", "C", "G", "T", "A", "C", "G", "T", "A", "C", "G", "T")[:n]
    ]


def test_build_array_1_has_1_fragment_no_stitcher() -> None:
    arr = build_array(crrnas=_dummy_crrnas(1))
    assert arr.array_size == 1
    assert len(arr.fragments) == 1
    assert arr.stitcher is None
    # Single fragment spans backbone directly: left=B5, right=B3.
    assert arr.fragments[0].left_overhang == constants.B5_OVERHANG
    assert arr.fragments[0].right_overhang == constants.B3_OVERHANG
    assert arr.total_pieces() == 1


def test_build_array_2_has_2_fragments_no_stitcher() -> None:
    arr = build_array(crrnas=_dummy_crrnas(2))
    assert arr.array_size == 2
    assert len(arr.fragments) == 2
    assert arr.stitcher is None
    # Position 1: B5 → J1; Position 2: J1 → B3.
    assert arr.fragments[0].left_overhang == constants.B5_OVERHANG
    assert arr.fragments[0].right_overhang == constants.J_OVERHANGS[0]
    assert arr.fragments[1].left_overhang == constants.J_OVERHANGS[0]
    assert arr.fragments[1].right_overhang == constants.B3_OVERHANG
    assert arr.total_pieces() == 2


def test_build_array_4_has_4_fragments_plus_stitcher() -> None:
    arr = build_array(crrnas=_dummy_crrnas(4))
    assert arr.array_size == 4
    assert len(arr.fragments) == 4
    assert arr.stitcher is not None
    assert arr.stitcher.stop_size == 4
    assert arr.total_pieces() == 5


def test_build_array_12_has_no_stitcher() -> None:
    arr = build_array(crrnas=_dummy_crrnas(12))
    assert arr.array_size == 12
    assert len(arr.fragments) == 12
    assert arr.stitcher is None
    assert arr.total_pieces() == 12


def test_build_array_fragment_overhangs_chain_correctly() -> None:
    """Adjacent fragments share the same overhang at their shared junction."""
    arr = build_array(crrnas=_dummy_crrnas(6))
    for i in range(len(arr.fragments) - 1):
        assert arr.fragments[i].right_overhang == arr.fragments[i + 1].left_overhang


def test_build_array_terminal_fragment_uses_stitcher_overhang_for_sub12() -> None:
    arr = build_array(crrnas=_dummy_crrnas(8))
    assert arr.stitcher is not None
    assert arr.fragments[-1].right_overhang == arr.stitcher.left_overhang
    assert arr.stitcher.right_overhang == constants.B3_OVERHANG


def test_build_array_terminal_fragment_uses_b3_for_size_12() -> None:
    arr = build_array(crrnas=_dummy_crrnas(12))
    assert arr.fragments[-1].right_overhang == constants.B3_OVERHANG


def test_build_array_barcode_only_on_position_1() -> None:
    arr = build_array(crrnas=_dummy_crrnas(4), barcode="CATGACGTACGT")
    assert arr.fragments[0].barcode == "CATGACGTACGT"
    for frag in arr.fragments[1:]:
        assert frag.barcode is None


def test_build_array_rejects_mismatched_crrna_count() -> None:
    with pytest.raises(ValueError, match="must equal array_size"):
        build_array(crrnas=_dummy_crrnas(3), array_size=4)


def test_build_stitcher_returns_none_for_sizes_with_b3_terminus() -> None:
    # Sizes 1, 2, 12 all have their terminal fragment go straight to B3 — no stitcher.
    assert build_stitcher(1) is None
    assert build_stitcher(2) is None
    assert build_stitcher(12) is None


def test_build_stitcher_overhangs_match_dossier() -> None:
    s = build_stitcher(4)
    assert s is not None
    assert s.left_overhang == constants.STITCHER_JUNCTION_BY_STOP_SIZE[4]
    assert s.right_overhang == constants.B3_OVERHANG
    assert constants.STITCHER_SPACER in s.ordered_dna
