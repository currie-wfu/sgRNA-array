"""Per-fragment and per-array assembly.

Generates the full ordered DNA sequences for each sgRNA fragment in an array and the
stitcher oligo for sub-12 array sizes. Returns structured `Fragment` / `Array` objects
that exporters consume.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from sgrna_array.constants import (
    B3_OVERHANG,
    B5_OVERHANG,
    J_OVERHANGS,
    SGRNA_SCAFFOLD_DEWEIRDT,
    SGRNA_SCAFFOLD_DEWEIRDT_CS1,
    STITCHER_JUNCTION_BY_STOP_SIZE,
    STITCHER_SPACER,
    SUPPORTED_ARRAY_SIZES,
    TERMINAL_USES_B3_DIRECTLY,
)
from sgrna_array.enzymes import wrap_paqci_cassette
from sgrna_array.ribozyme import build_full_hh, build_hdv


@dataclass(frozen=True)
class Fragment:
    """A single ordered DNA fragment encoding one ribozyme-flanked sgRNA."""

    position: int  # 1-based position in the array
    gene_label: str  # user-supplied label (e.g., target gene name)
    crrna: str  # 20-nt spacer
    left_overhang: str  # 4-nt 5' overhang revealed after PaqCI digestion
    right_overhang: str  # 4-nt 3' overhang
    barcode: str | None  # 5' barcode (only on position-1 fragment)
    use_cs1: bool  # whether the scaffold has CS1 in Hairpin-1 loop
    ordered_dna: str  # the full DNA sequence to order

    def __len__(self) -> int:
        return len(self.ordered_dna)


@dataclass(frozen=True)
class Stitcher:
    """A short stitcher oligo that terminates an array before position 12."""

    stop_size: int  # which sub-size this stitcher terminates (4, 6, 8, or 10)
    left_overhang: str  # matches Jₙ where n = stop_size
    right_overhang: str  # = B3
    ordered_dna: str


@dataclass(frozen=True)
class Array:
    """The full ordered set of pieces (Fragments + optional Stitcher) for one array design."""

    array_size: int  # 4, 6, 8, 10, or 12
    fragments: list[Fragment] = field(default_factory=list)
    stitcher: Stitcher | None = None

    def total_pieces(self) -> int:
        """Number of separate DNA pieces to order for this array."""
        return len(self.fragments) + (1 if self.stitcher is not None else 0)


def junction_overhang(position: int, array_size: int) -> str:
    """Return the 3' overhang for the fragment at this position.

    The 3' overhang of fragment N is the same physical 4-nt sequence as the 5' overhang
    of fragment N+1 (they form the ligation junction). For the terminal position of the
    array, the 3' overhang is B3 directly (sizes 1, 2, 12 — see
    `TERMINAL_USES_B3_DIRECTLY`) or the Jₙ overhang that a stitcher oligo bridges to B3
    (sizes 4, 6, 8, 10).
    """
    if position == array_size:
        if array_size in TERMINAL_USES_B3_DIRECTLY:
            return B3_OVERHANG
        return STITCHER_JUNCTION_BY_STOP_SIZE[array_size]
    # Non-terminal: 3' overhang is the junction to position + 1.
    return J_OVERHANGS[position - 1]  # J1 is at index 0


def left_overhang(position: int) -> str:
    """Return the 5' overhang for the fragment at this position.

    Position 1's 5' overhang is B5 (backbone). Position N>1's 5' overhang is the same as
    position N-1's 3' overhang (the shared junction).
    """
    if position == 1:
        return B5_OVERHANG
    return J_OVERHANGS[position - 2]  # J(N-1) is at index N-2


def build_fragment(
    *,
    crrna: str,
    position: int,
    array_size: int,
    gene_label: str = "",
    barcode: str | None = None,
    use_cs1: bool = False,
) -> Fragment:
    """Build a single Fragment for the given crRNA and array position.

    Args:
        crrna: 20-nt spacer sequence. Should already have passed `validator.validate_crrna`.
        position: 1-based position in the array (1 .. array_size).
        array_size: One of SUPPORTED_ARRAY_SIZES.
        gene_label: Optional human-readable label carried through to the export.
        barcode: 5' array barcode (only attached if `position == 1`). If None and
            position == 1, no barcode is prepended. Callers normally generate or
            validate the barcode via `barcode` module before passing it here.
        use_cs1: If True, use the CS1-modified DeWeirdt scaffold (10x direct capture);
            otherwise use the default DeWeirdt scaffold.

    Returns:
        A `Fragment` whose `ordered_dna` is the full DNA sequence ready to send to the
        oligo / gBlock vendor.
    """
    _validate_position(position, array_size)
    scaffold = SGRNA_SCAFFOLD_DEWEIRDT_CS1 if use_cs1 else SGRNA_SCAFFOLD_DEWEIRDT
    crrna_lc = crrna.lower()

    # Build the RNA-coding core (stem-I + HH core + crRNA + scaffold + HDV).
    core_segments = [build_full_hh(crrna_lc), crrna_lc, scaffold, build_hdv()]

    # Position-1 fragment gets the 5' barcode prepended to the core.
    if position == 1 and barcode:
        core_segments.insert(0, barcode.upper())

    core = "".join(core_segments)

    left = left_overhang(position)
    right = junction_overhang(position, array_size)

    ordered = wrap_paqci_cassette(core=core, left_overhang=left, right_overhang=right)

    return Fragment(
        position=position,
        gene_label=gene_label,
        crrna=crrna_lc,
        left_overhang=left,
        right_overhang=right,
        barcode=barcode if position == 1 else None,
        use_cs1=use_cs1,
        ordered_dna=ordered,
    )


def build_stitcher(array_size: int) -> Stitcher | None:
    """Build the stitcher oligo for sizes that need one; None for sizes 1, 2, 12.

    Sizes 1, 2, and 12 have their terminal fragment ligate to backbone B3 directly, so
    no stitcher is required. Sizes 4, 6, 8, 10 use a library-style stitcher that bridges
    the Jₙ junction to B3.
    """
    if array_size in TERMINAL_USES_B3_DIRECTLY:
        return None
    if array_size not in STITCHER_JUNCTION_BY_STOP_SIZE:
        raise ValueError(
            f"No stitcher defined for array_size={array_size}; "
            f"sizes needing a stitcher: {sorted(STITCHER_JUNCTION_BY_STOP_SIZE)}"
        )
    left = STITCHER_JUNCTION_BY_STOP_SIZE[array_size]
    right = B3_OVERHANG
    # Stitchers are short and don't carry PaqCI sites; they are short duplex oligos with
    # the matching 5' overhangs ordered as ssDNA pairs that anneal during the reaction.
    # The "ordered DNA" representation here is the top strand of the duplex region.
    ordered = left + STITCHER_SPACER + right
    return Stitcher(
        stop_size=array_size,
        left_overhang=left,
        right_overhang=right,
        ordered_dna=ordered,
    )


def build_array(
    *,
    crrnas: list[str],
    gene_labels: list[str] | None = None,
    array_size: int | None = None,
    barcode: str | None = None,
    use_cs1: bool = False,
) -> Array:
    """Build a full Array from a list of crRNAs (one per position).

    Args:
        crrnas: List of 20-nt spacers, ordered by array position (index 0 = position 1).
        gene_labels: Optional parallel list of gene labels; defaults to empty strings.
        array_size: Defaults to the length of `crrnas`; must be in SUPPORTED_ARRAY_SIZES
            and ≥ len(crrnas).
        barcode: Optional 5' barcode attached to the position-1 fragment.
        use_cs1: If True, all fragments use the CS1-modified scaffold.

    Returns:
        An `Array` with `len(crrnas)` Fragments plus an optional Stitcher for sub-12 sizes.
    """
    if not crrnas:
        raise ValueError("At least one crRNA required")
    array_size = array_size or len(crrnas)
    if array_size not in SUPPORTED_ARRAY_SIZES:
        raise ValueError(
            f"array_size={array_size} not in supported sizes {SUPPORTED_ARRAY_SIZES}"
        )
    if len(crrnas) != array_size:
        raise ValueError(
            f"Number of crRNAs ({len(crrnas)}) must equal array_size ({array_size})"
        )
    labels = gene_labels or [""] * array_size
    if len(labels) != array_size:
        raise ValueError(
            f"gene_labels length ({len(labels)}) must equal array_size ({array_size})"
        )

    fragments = [
        build_fragment(
            crrna=crrna,
            position=i + 1,
            array_size=array_size,
            gene_label=labels[i],
            barcode=barcode if i == 0 else None,
            use_cs1=use_cs1,
        )
        for i, crrna in enumerate(crrnas)
    ]

    stitcher = build_stitcher(array_size)

    return Array(array_size=array_size, fragments=fragments, stitcher=stitcher)


def _validate_position(position: int, array_size: int) -> None:
    if array_size not in SUPPORTED_ARRAY_SIZES:
        raise ValueError(
            f"array_size={array_size} not in supported sizes {SUPPORTED_ARRAY_SIZES}"
        )
    if not (1 <= position <= array_size):
        raise ValueError(
            f"position={position} out of range for array_size={array_size}"
        )
