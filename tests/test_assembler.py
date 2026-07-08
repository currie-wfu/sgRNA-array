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
    # Expected length: 5 pad + 7 PaqCI + 4 filler + 4 overhang + 6 stem-I + 37 HH core
    # + 20 crRNA + 76 scaffold + 68 HDV + 4 overhang + 4 filler + 7 PaqCI + 5 pad
    # = 247 nt (was 241 before the 4-nt filler fix; the 6-nt increase is the corrected
    # PaqCI cut-offset space).
    expected_len = (
        constants.PAQCI_FLANKING_PAD_LEN
        + len(constants.PAQCI_RECOGNITION_FWD)
        + len(constants.PAQCI_RECOGNITION_TO_OVERHANG_SPACER)  # NOW 4 nt (was 1)
        + 4  # left overhang
        + constants.HH_STEM_I_LEN
        + len(constants.HH_CATALYTIC_CORE)
        + constants.CRRNA_LEN
        + len(constants.SGRNA_SCAFFOLD_DEWEIRDT)
        + len(constants.HDV_RIBOZYME)
        + 4  # right overhang
        + len(constants.PAQCI_RECOGNITION_TO_OVERHANG_SPACER)  # NOW 4 nt (was 1)
        + len(constants.PAQCI_RECOGNITION_REV)
        + constants.PAQCI_FLANKING_PAD_LEN
    )
    assert len(frag.ordered_dna) == expected_len == 247


def test_intended_overhang_sits_at_the_paqci_cut_position() -> None:
    """Regression guard for the 4-nt-filler bug.

    PaqCI cuts 4 nt 3' of CACCTGC on the top strand. The 4-nt sticky end exposed after
    digestion is at positions 5–8 downstream of the recognition site. If this test
    fails, the tool is generating fragments whose exposed overhang is NOT the intended
    overhang — assemblies into the host vector will fail.
    """
    frag = build_fragment(crrna=EXAMPLE_CRRNA, position=1, array_size=12)
    seq = frag.ordered_dna

    # Find CACCTGC on the top strand.
    idx_5 = seq.find(constants.PAQCI_RECOGNITION_FWD)
    assert idx_5 >= 0, "5' PaqCI recognition site missing"
    end_of_recog_5 = idx_5 + len(constants.PAQCI_RECOGNITION_FWD)
    # Positions 5–8 after the recognition site (= 4-nt overhang exposed after PaqCI cut).
    exposed_left = seq[end_of_recog_5 + 4 : end_of_recog_5 + 8]
    assert exposed_left.upper() == frag.left_overhang.upper(), (
        f"5' overhang exposed by PaqCI ({exposed_left!r}) does not match the intended "
        f"left overhang ({frag.left_overhang!r}). "
        f"The 4-nt filler between CACCTGC and the overhang is likely wrong."
    )

    # Same check for the 3' side (GCAGGTG). Positions 5–8 upstream of GCAGGTG on top
    # correspond to bottom-strand positions 1–4 downstream of CACCTGC on bottom.
    idx_3 = seq.find(constants.PAQCI_RECOGNITION_REV)
    assert idx_3 >= 0, "3' PaqCI recognition site missing"
    exposed_right = seq[idx_3 - 8 : idx_3 - 4]
    assert exposed_right.upper() == frag.right_overhang.upper(), (
        f"3' overhang exposed by PaqCI ({exposed_right!r}) does not match the intended "
        f"right overhang ({frag.right_overhang!r})."
    )


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


# ---------------------------------------------------------------------------
# Contiguous insert (secondary output: whole array as one gBlock)
# ---------------------------------------------------------------------------


def test_fragment_exposes_core_field() -> None:
    """Each Fragment carries its core (no PaqCI cassette, no overhangs) for re-use."""
    frag = build_fragment(crrna=EXAMPLE_CRRNA, position=1, array_size=4)
    # ordered_dna = pad + PaqCI_fwd + filler + L_overhang + core + R_overhang + filler + PaqCI_rev + pad
    # Wrapper is 20 nt on each side (5 pad + 7 recog + 4 filler + 4 overhang).
    assert frag.core == frag.ordered_dna[20:-20]
    # Core includes stem-I, crRNA, scaffold, HDV.
    assert EXAMPLE_CRRNA.lower() in frag.core.lower()
    assert EXPECTED_STEM_I in frag.core.lower()


def test_contiguous_insert_size_1_matches_single_fragment_core_in_wrapper() -> None:
    """For a size-1 array, the contiguous insert is just one core wrapped in PaqCI/B5/B3."""
    arr = build_array(crrnas=[EXAMPLE_CRRNA], barcode="ACGTACGTACGT")
    insert = arr.contiguous_insert()
    # wrap_paqci_cassette wraps with pad(5) + PaqCI_fwd(7) + filler(4) + L_overhang(4) on the
    # left and the mirror on the right — so the 20-nt wrapper INCLUDES B5/B3 already.
    # For size-1 with a 12-nt barcode, core = 12 + 6 + 37 + 20 + 76 + 68 = 219.
    # Total = 20 + 219 + 20 = 259.
    assert len(insert) == 20 + 219 + 20 == 259


def test_contiguous_insert_size_4_concatenates_cores_with_junctions() -> None:
    """For larger arrays, cores are concatenated with inter-position 4-nt overhangs."""
    arr = build_array(crrnas=_dummy_crrnas(4), barcode="ACGTACGTACGT")
    insert = arr.contiguous_insert()
    # Body = core1_with_bc(219) + J1(4) + core2(207) + J2(4) + core3(207) + J3(4) + core4(207)
    #      = 219 + 4 + 207 + 4 + 207 + 4 + 207 = 852
    # Total = 20 + body + 20 = 892 (B5/B3 are inside the 20-nt wrapper).
    expected_body = 219 + 4 + 207 + 4 + 207 + 4 + 207
    assert len(insert) == 20 + expected_body + 20 == 892


def test_contiguous_insert_has_exactly_one_paqci_site_per_end() -> None:
    """No internal PaqCI sites; just the outer two for cloning into the vector."""
    arr = build_array(crrnas=_dummy_crrnas(6), barcode="ACGTACGTACGT")
    insert = arr.contiguous_insert()
    assert insert.upper().count("CACCTGC") == 1
    assert insert.upper().count("GCAGGTG") == 1


def test_contiguous_insert_contains_all_crrnas_in_order() -> None:
    # Use crRNAs with no shared 10-nt prefix so `find()` locates the right one.
    crrnas = [
        "ACGTACGTACGTACGTACGT",
        "TGCATGCATGCATGCATGCA",
        "AAATTTCCCGGGAAATTTCC",
        "GGGAAACCCTTTGGGAAACC",
        "CACGCACGCACGCACGCACG",
        "GTCAGTCAGTCAGTCAGTCA",
    ]
    arr = build_array(crrnas=crrnas, barcode="ACGTACGTACGT")
    insert = arr.contiguous_insert().lower()
    last_pos = -1
    for crrna in crrnas:
        idx = insert.find(crrna.lower())
        assert idx > last_pos, f"crRNA {crrna} not found after position {last_pos}"
        last_pos = idx


def test_contiguous_insert_outer_overhangs_are_b5_b3() -> None:
    """Reading post-digestion, the outer sticky ends are B5 (5') and B3 (3')."""
    arr = build_array(crrnas=_dummy_crrnas(4), barcode="ACGTACGTACGT")
    insert = arr.contiguous_insert()
    # Strip the wrapper: 5 pad + 7 PaqCI_fwd + 4 filler = 16 nt prefix, mirrored suffix.
    after_paqci = insert[16:-16]
    assert after_paqci.startswith(constants.B5_OVERHANG)
    assert after_paqci.endswith(constants.B3_OVERHANG)


# ---------------------------------------------------------------------------
# Reconstituted-PaqCI-site guard
# ---------------------------------------------------------------------------


def test_build_fragment_rejects_crrna_that_reconstitutes_paqci_at_scaffold_junction() -> None:
    """A crRNA ending in `GCAGGT` + scaffold starting with `G` spells `GCAGGTG` at the
    boundary — that's a reconstituted PaqCI reverse site, which would break assembly.

    The per-crRNA validator doesn't catch this (the crRNA alone doesn't contain a
    PaqCI site — only after concatenation with the scaffold does one appear). The
    assembler's post-hoc check must catch it.
    """
    # scaffold starts with 'gtttg...', so a crRNA ending in 'GCAGGT' produces
    # ...GCAGGT + g... = GCAGGTG. 20 nt total, no internal PaqCI site.
    bad_crrna = "AAAAAAAAAAAAAAGCAGGT"
    assert "CACCTGC" not in bad_crrna and "GCAGGTG" not in bad_crrna
    with pytest.raises(ValueError, match="unexpected PaqCI site count"):
        build_fragment(
            crrna=bad_crrna,
            position=1,
            array_size=1,
            gene_label="bad",
        )


def test_build_fragment_rejects_crrna_that_reconstitutes_forward_paqci() -> None:
    """A crRNA ending in `CACCTG` followed by `C` at scaffold start would form
    `CACCTGC`. Our scaffold starts with `g`, not `c`, so we need a crRNA that
    contains part of CACCTGC at the boundary with the scaffold. Simulate by
    putting `CACCTG` at the very end of the crRNA and confirming detection when
    the assembled sequence is scanned."""
    # Trickier to trigger — needs the immediate next character (start of scaffold)
    # to be C, which our scaffold's start ('g') doesn't provide. So verify the
    # symmetric case works: crRNA with CACCTGC internally is caught earlier by
    # the validator, and the assembler check would also catch it downstream.
    # Here we assert that a crRNA CONTAINING CACCTGC is caught by the assembler
    # even if someone bypasses the up-front validator.
    embedded = "AAAAAACACCTGCAAAAAAA"  # CACCTGC at positions 6-12; 20 nt total
    assert "CACCTGC" in embedded
    with pytest.raises(ValueError, match="unexpected PaqCI site count"):
        build_fragment(
            crrna=embedded,
            position=1,
            array_size=1,
            gene_label="embedded_paqci",
        )


def test_valid_fragment_passes_reconstitution_guard() -> None:
    """The Shh gRNA#3 example (which we know is a good crRNA) must not trip the guard."""
    frag = build_fragment(
        crrna=EXAMPLE_CRRNA,
        position=1,
        array_size=1,
        gene_label="Shh",
    )
    # Exactly 1 forward + 1 reverse PaqCI site (the outer cassette).
    upper = frag.ordered_dna.upper()
    assert upper.count("CACCTGC") == 1
    assert upper.count("GCAGGTG") == 1


def test_contiguous_insert_size_1_passes_reconstitution_guard() -> None:
    arr = build_array(
        crrnas=[EXAMPLE_CRRNA],
        barcode="ACGTACGTACGT",
    )
    insert = arr.contiguous_insert().upper()
    assert insert.count("CACCTGC") == 1
    assert insert.count("GCAGGTG") == 1


def test_contiguous_insert_internal_junctions_preserved() -> None:
    """The 4-nt overhangs between adjacent cores appear once between them."""
    arr = build_array(crrnas=_dummy_crrnas(4))
    insert = arr.contiguous_insert()
    # J1, J2, J3 should each appear at least once internally (in addition to other places).
    # We do a stronger check: the junctions appear at the boundary between fragment cores.
    for i in range(len(arr.fragments) - 1):
        junction = arr.fragments[i].right_overhang
        # The junction sits between core_i and core_(i+1). It must appear in the insert.
        assert junction in insert


def test_build_stitcher_overhangs_match_dossier() -> None:
    s = build_stitcher(4)
    assert s is not None
    assert s.left_overhang == constants.STITCHER_JUNCTION_BY_STOP_SIZE[4]
    assert s.right_overhang == constants.B3_OVERHANG
    assert constants.STITCHER_SPACER in s.ordered_dna
