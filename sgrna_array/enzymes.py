"""PaqCI Type IIs enzyme handling: site scanning and cassette construction.

PaqCI recognition: `CACCTGC` (forward) / `GCAGGTG` (reverse).
Cuts 4 nt 3' of recognition on the top strand, 8 nt 3' on the bottom strand,
leaving a 4-nt 5' overhang.

The per-fragment cassette wraps the [stem-I + HH + crRNA + scaffold + HDV] RNA-coding
core with PaqCI recognition sites pointing inward, so that one-pot Golden Gate
digestion releases the core with the position-specific 4-nt overhangs revealed.
"""

from __future__ import annotations

from sgrna_array.constants import (
    PAQCI_FLANKING_PAD_LEN,
    PAQCI_RECOGNITION_FWD,
    PAQCI_RECOGNITION_REV,
    PAQCI_RECOGNITION_TO_OVERHANG_SPACER,
)

#: Default neutral padding sequence to flank the PaqCI cassette. Length matches
#: PAQCI_FLANKING_PAD_LEN; chosen to be inert (no splice motifs, no enzyme sites).
_DEFAULT_PAD: str = "ATACT"  # 5 nt
assert len(_DEFAULT_PAD) == PAQCI_FLANKING_PAD_LEN


def scan_paqci_sites(seq: str) -> list[tuple[int, str]]:
    """Find PaqCI recognition sites in a DNA sequence (both orientations).

    Returns a list of `(start_position, orientation)` tuples, where orientation is
    `'fwd'` or `'rev'`. Used by `validator` to reject crRNAs containing internal
    PaqCI sites that would break Golden Gate assembly.
    """
    seq_upper = seq.upper()
    hits: list[tuple[int, str]] = []
    for i in range(len(seq_upper) - len(PAQCI_RECOGNITION_FWD) + 1):
        window = seq_upper[i : i + len(PAQCI_RECOGNITION_FWD)]
        if window == PAQCI_RECOGNITION_FWD:
            hits.append((i, "fwd"))
        elif window == PAQCI_RECOGNITION_REV:
            hits.append((i, "rev"))
    return hits


def wrap_paqci_cassette(
    core: str,
    left_overhang: str,
    right_overhang: str,
    pad: str = _DEFAULT_PAD,
) -> str:
    """Wrap an RNA-coding core in a PaqCI Golden Gate cassette.

    Layout (5'→3', top strand):
        [pad] - CACCTGC - N - [left_overhang] - [core] - [right_overhang] - N - GCAGGTG - [pad]

    When PaqCI digests this fragment, the recognition sites + flanking N + pad are removed,
    releasing the core with the two 4-nt overhangs exposed as 5' sticky ends.

    Args:
        core: The RNA-coding sequence (stem-I + HH + crRNA + scaffold + HDV).
        left_overhang: 4-nt overhang to expose at the 5' end of the core after digestion.
        right_overhang: 4-nt overhang at the 3' end.
        pad: Neutral flanking sequence outside each recognition site. Defaults to
            `_DEFAULT_PAD`.

    Returns:
        The full ordered DNA sequence for this fragment.
    """
    if len(left_overhang) != 4 or len(right_overhang) != 4:
        raise ValueError(
            f"Overhangs must be 4 nt; got left={len(left_overhang)}, right={len(right_overhang)}"
        )
    spacer = PAQCI_RECOGNITION_TO_OVERHANG_SPACER
    return (
        pad
        + PAQCI_RECOGNITION_FWD
        + spacer
        + left_overhang
        + core
        + right_overhang
        + spacer
        + PAQCI_RECOGNITION_REV
        + pad
    )
