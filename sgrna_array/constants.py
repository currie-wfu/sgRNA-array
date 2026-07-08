"""Locked sequence constants for the sgRNA array webtool.

Every sequence here is verified against a cited source — either a peer-reviewed
publication, a deposited Addgene plasmid, or the host vector GenBank file.
DO NOT modify these without updating the dossier in memory and the citations below.

Key references:
- DeWeirdt et al. 2022, Nat Commun, DOI 10.1038/s41467-022-33024-2 (sgRNA scaffold)
- Replogle et al. 2020, Nat Biotechnol, DOI 10.1038/s41587-020-0470-y (CS1, CS2)
- 10x Genomics CG000197 RevA (CS1, CS2 cassette positions)
- Potapov et al. 2018, ACS Synth Biol (T4 ligase 4-nt overhang fidelity)
- NEB Golden Gate Assembly Tool generate-mode at goldengate.neb.com (overhang set, 98% fidelity)
- Host vector: 1103p_SceITol2_14xUAS_CRISPR2-0ShhmNeonGreen-T2A-Cas9.gb
- pRDA_118 / Addgene #133459 (DeWeirdt scaffold deposit)
"""

from __future__ import annotations

# ---------------------------------------------------------------------------
# Type IIs enzyme
# ---------------------------------------------------------------------------

#: PaqCI recognition site on the top strand (5'→3').
PAQCI_RECOGNITION_FWD: str = "CACCTGC"

#: PaqCI recognition site on the top strand when the recognition is on the bottom (reverse).
PAQCI_RECOGNITION_REV: str = "GCAGGTG"

#: PaqCI cut offsets relative to the 3' end of the recognition site. Top-strand cut is
#: 4 nt 3' of the recognition; bottom-strand cut is 8 nt 3'. Net result: 4-nt 5' overhang.
PAQCI_CUT_OFFSET_TOP: int = 4
PAQCI_CUT_OFFSET_BOTTOM: int = 8

#: Number of "padding" nucleotides outside the PaqCI recognition site on each end of an
#: ordered DNA fragment. Sets how much single-stranded DNA sits between the last
#: PaqCI cut and the very end of the ordered piece.
#:
#: NEB's Golden Gate FAQ for PaqCI (2024-12-10 revision) explicitly recommends
#: **6 flanking bases** at the 5' ends of Golden Gate amplicons: "For Golden Gate
#: Assembly, how many base pairs should my amplicon inserts have flanking the PaqCI
#: restriction site?" — answer: 6. This is above PaqCI's general cleavage-near-the-end
#: profile (50–100% activity at ≥2 nt from a DNA terminus; 20–50% at 1 nt; 0% at 0).
#: One-pot Golden Gate needs maximum enzyme activity to drive assembly to completion,
#: hence the ≥6 nt guidance (vs. the ≥2 nt hard floor).
#:
#: Bumping from 5 → 6 nt costs +2 nt per fragment (+2 per contiguous insert).
PAQCI_FLANKING_PAD_LEN: int = 6

#: 4-nt filler between the PaqCI recognition site and the 4-nt overhang.
#:
#: PaqCI is CACCTGC(4/8) — top-strand cut 4 nt 3' of recognition, bottom-strand cut 8 nt.
#: So the actual 4-nt overhang exposed after digestion sits at positions 5–8 downstream
#: of the recognition site. Positions 1–4 (this filler) get cut away with the recognition
#: site and are discarded.
#:
#: WITHOUT this filler (or with fewer than 4 nt of it), PaqCI cuts INSIDE the intended
#: overhang and the exposed 4-nt sticky end is (nt 4 of intended overhang) + (first 3 nt
#: of whatever follows). That's a bug — the exposed overhang no longer matches what the
#: backbone expects, and downstream contents (barcode, stem-I) leak into it and change
#: from design to design.
#:
#: Any 4 nt works — this sequence gets discarded during digestion. `AATA` chosen for
#: 25 % GC, non-palindromic, and it doesn't reconstitute a PaqCI site when concatenated
#: with `CACCTGC` or `GCAGGTG`.
PAQCI_RECOGNITION_TO_OVERHANG_SPACER: str = "AATA"

# Sanity check: this filler MUST match the top-strand cut offset, or PaqCI cuts in the
# wrong place. If someone edits the constant above, the assertion below will fire at
# module import time and prevent a bad build from ever running.
assert len(PAQCI_RECOGNITION_TO_OVERHANG_SPACER) == PAQCI_CUT_OFFSET_TOP, (
    f"PAQCI_RECOGNITION_TO_OVERHANG_SPACER must be {PAQCI_CUT_OFFSET_TOP} nt to place "
    f"the intended overhang at the correct PaqCI cut site; got "
    f"{len(PAQCI_RECOGNITION_TO_OVERHANG_SPACER)} nt "
    f"({PAQCI_RECOGNITION_TO_OVERHANG_SPACER!r})"
)


# ---------------------------------------------------------------------------
# Hammerhead (HH) ribozyme
# ---------------------------------------------------------------------------

#: HH catalytic core (37 nt, constant). Sits between the per-crRNA stem-I and the crRNA spacer
#: in the encoded fragment. Verified against host vector positions 1928-1964.
HH_CATALYTIC_CORE: str = "ctgatgagtccgtgaggacgaaacgagtaagctcgtc"

#: Length of the variable stem-I region (per crRNA). Stem-I = reverse complement of crRNA[:6],
#: forming a 6-bp Stem I helix with the spacer.
HH_STEM_I_LEN: int = 6


# ---------------------------------------------------------------------------
# HDV ribozyme
# ---------------------------------------------------------------------------

#: HDV ribozyme (68 nt, constant). Sits 3' of the sgRNA scaffold and self-cleaves to give a
#: defined 3' end on the released mature sgRNA. Verified against host vector positions 2065-2132.
HDV_RIBOZYME: str = (
    "ggccggcatggtcccagcctcctcgctggcgccggctgggcaacatgcttcggcatggcgaatgggac"
)


# ---------------------------------------------------------------------------
# sgRNA scaffold
# ---------------------------------------------------------------------------

#: DeWeirdt 2022 sgRNA scaffold (76 nt, no Pol III terminator).
#: Verified from pRDA_118 / Addgene #133459. Hsu/Mali backbone with T5→G and A26→C
#: (compensatory mutation restoring the lower-stem 5↔26 base pair as G:C).
#: Used as the default scaffold for new arrays.
SGRNA_SCAFFOLD_DEWEIRDT: str = (
    "gtttgagagctagaaatagcaagttcaaataaggctagtccgttatcaacttgaaaaagtggcaccgagtcggtgc"
)

#: CS1-modified DeWeirdt scaffold (102 nt). Replaces the 4-nt Hairpin-1 loop `gaaa`
#: (positions 53-56 of SGRNA_SCAFFOLD_DEWEIRDT) with `ggcc-CS1-ggcc` (4+22+4 = 30 nt).
#: Implements the Replogle 2020 `sgRNA-CR1cs1` modification on the DeWeirdt backbone.
SGRNA_SCAFFOLD_DEWEIRDT_CS1: str = (
    "gtttgagagctagaaatagcaagttcaaataaggctagtccgttatcaactt"
    "ggccgctttaaggccggtcctagcaaggcc"
    "aagtggcaccgagtcggtgc"
)

#: Hsu/Mali 2013 scaffold (76 nt). Retained for reference and vector-consistency mode;
#: not the v1 default.
SGRNA_SCAFFOLD_HSU_MALI: str = (
    "gttttagagctagaaatagcaagttaaaataaggctagtccgttatcaacttgaaaaagtggcaccgagtcggtgc"
)


# ---------------------------------------------------------------------------
# Illumina Capture Sequences (Replogle 2020 / 10x CG000197 RevA)
# ---------------------------------------------------------------------------

#: Capture Sequence 1, in-sgRNA orientation (22 nt). Annealing partner on the 10x gel bead
#: is the reverse complement (`TTGCTAGGACCGGCCTTAAAGC`). Must be placed in the Hairpin-1
#: loop; 3'-end placement compromises guide activity.
CS1_INSGRNA: str = "GCTTTAAGGCCGGTCCTAGCAA"

#: Capture Sequence 1, on-gel-bead orientation (22 nt; reverse complement of CS1_INSGRNA).
CS1_GELBEAD: str = "TTGCTAGGACCGGCCTTAAAGC"

#: Capture Sequence 2, in-sgRNA orientation (22 nt). Not exposed in v1; flexible placement.
CS2_INSGRNA: str = "GCTCACCTATTAGCGGCTAAGG"

#: Capture Sequence 2, on-gel-bead orientation (22 nt).
CS2_GELBEAD: str = "CCTTAGCCGCTAATAGGTGAGC"

#: The 4-nt closing stem flanking CS1/CS2 inside the modified Hairpin-1 loop.
#: Self-palindromic, base-pairs with itself.
CS_HAIRPIN_STEM: str = "ggcc"


# ---------------------------------------------------------------------------
# Golden Gate overhang set — NEB-validated at 98% fidelity
# ---------------------------------------------------------------------------

#: Backbone-to-position-1 overhang. Locked by the host vector's 5' PaqCI rev cut
#: (vector positions 1907-1910).
B5_OVERHANG: str = "ACGG"

#: Position-12-to-backbone overhang. Locked by the host vector's 3' PaqCI fwd cut
#: (vector positions 2144-2147).
B3_OVERHANG: str = "GAGC"

#: Inter-position junction overhangs J1…J11. Order matters: J[i] is the overhang at the
#: junction between position i and position i+1. Validated together with B5/B3 by NEB's
#: Golden Gate Assembly Tool generate-mode (98% predicted fidelity).
J_OVERHANGS: tuple[str, ...] = (
    "ATAA",  # J1: position 1 → 2
    "AAAT",  # J2: position 2 → 3
    "CTCA",  # J3: position 3 → 4
    "ACCA",  # J4: position 4 → 5  (also = stop-at-4 stitcher overhang)
    "CCAA",  # J5: position 5 → 6
    "CAGA",  # J6: position 6 → 7  (also = stop-at-6 stitcher overhang)
    "GAAA",  # J7: position 7 → 8
    "AGTA",  # J8: position 8 → 9  (also = stop-at-8 stitcher overhang)
    "AATG",  # J9: position 9 → 10
    "AGAC",  # J10: position 10 → 11 (also = stop-at-10 stitcher overhang)
    "CTAC",  # J11: position 11 → 12
)

#: Map from "stop-at-N" array size to the junction overhang whose stitcher oligo terminates
#: the assembly. The stitcher is a short duplex with the J_n overhang on one end and B3 on
#: the other, terminating the array at position N.
STITCHER_JUNCTION_BY_STOP_SIZE: dict[int, str] = {
    4: J_OVERHANGS[3],   # J4 = ACCA
    6: J_OVERHANGS[5],   # J6 = CAGA
    8: J_OVERHANGS[7],   # J8 = AGTA
    10: J_OVERHANGS[9],  # J10 = AGAC
    # 12: natural endpoint at B3; no stitcher needed.
    # 1, 2: terminal fragment's 3' overhang is B3 directly; no stitcher needed
    # (see TERMINAL_USES_B3_DIRECTLY below).
}

#: Inert spacer sequence between the two overhangs of a stitcher oligo. 6 nt of neutral
#: sequence that contains no PaqCI sites and no obvious splice motifs.
STITCHER_SPACER: str = "ATAATA"

#: Array sizes whose terminal-position fragment's 3' overhang is B3 directly, with no
#: stitcher needed. Size 12 is the natural endpoint of the full library. Sizes 1 and 2
#: are small enough that the terminal fragment can be ordered as a self-contained piece
#: that ligates straight into the backbone — no library-stitcher dependency.
TERMINAL_USES_B3_DIRECTLY: frozenset[int] = frozenset({1, 2, 12})

#: All valid array sizes the tool supports.
SUPPORTED_ARRAY_SIZES: tuple[int, ...] = (1, 2, 4, 6, 8, 10, 12)


# ---------------------------------------------------------------------------
# 5' array barcode
# ---------------------------------------------------------------------------

#: Default length of the per-array 5' barcode that identifies the array in bulk PCR/Sanger.
DEFAULT_BARCODE_LEN: int = 12

#: Minimum pairwise Hamming distance between barcodes generated within a single batch.
DEFAULT_BARCODE_MIN_HAMMING: int = 3

#: GC-content bounds for generated barcodes (inclusive). Range 0.0–1.0.
DEFAULT_BARCODE_GC_MIN: float = 0.40
DEFAULT_BARCODE_GC_MAX: float = 0.60

#: Maximum allowed homopolymer run in a barcode.
DEFAULT_BARCODE_MAX_HOMOPOLYMER: int = 2


# ---------------------------------------------------------------------------
# crRNA validation parameters
# ---------------------------------------------------------------------------

#: Required length of an input crRNA spacer in nt.
CRRNA_LEN: int = 20

#: Allowed alphabet for crRNAs and barcodes.
DNA_ALPHABET: frozenset[str] = frozenset("ACGTacgt")

#: GC-content soft-warn bounds for crRNAs (inclusive). Range 0.0–1.0.
CRRNA_GC_WARN_MIN: float = 0.25
CRRNA_GC_WARN_MAX: float = 0.75

#: Maximum homopolymer run allowed in a crRNA before a soft warning fires.
CRRNA_HOMOPOLYMER_WARN_THRESHOLD: int = 5


# ---------------------------------------------------------------------------
# Host vector reference (for documentation / validation only)
# ---------------------------------------------------------------------------

#: Host vector filename for reference. Tool does not require this file at runtime; it is
#: used to document what backbone the locked B5/B3 overhangs come from.
HOST_VECTOR_FILENAME: str = "1103p_SceITol2_14xUAS_CRISPR2-0ShhmNeonGreen-T2A-Cas9.gb"

#: Vector positions of the locked PaqCI sites and the resulting backbone overhangs.
HOST_VECTOR_PAQCI_5_POSITION: tuple[int, int] = (1915, 1921)  # GCAGGTG
HOST_VECTOR_PAQCI_3_POSITION: tuple[int, int] = (2133, 2139)  # CACCTGC
HOST_VECTOR_B5_OVERHANG_POSITION: tuple[int, int] = (1907, 1910)
HOST_VECTOR_B3_OVERHANG_POSITION: tuple[int, int] = (2144, 2147)
