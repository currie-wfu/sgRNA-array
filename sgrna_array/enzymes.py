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
#: PAQCI_FLANKING_PAD_LEN (6 nt per NEB's Golden Gate FAQ). Chosen to be inert:
#: no splice donor/acceptor motifs, no BsaI/BsmBI/PaqCI recognition when
#: concatenated with either PaqCI recognition site or with any of our overhangs.
_DEFAULT_PAD: str = "AATACT"  # 6 nt
assert len(_DEFAULT_PAD) == PAQCI_FLANKING_PAD_LEN, (
    f"_DEFAULT_PAD length ({len(_DEFAULT_PAD)}) must match "
    f"PAQCI_FLANKING_PAD_LEN ({PAQCI_FLANKING_PAD_LEN}). Update both in sync."
)


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


def assert_no_extra_paqci_sites(
    sequence: str,
    *,
    expected_fwd: int = 1,
    expected_rev: int = 1,
    context: str = "",
) -> None:
    """Raise ValueError if the sequence has an unexpected number of PaqCI sites.

    A well-formed ordered DNA fragment (or contiguous insert) has exactly one
    forward PaqCI site (`CACCTGC`) at the 5' cassette and one reverse site
    (`GCAGGTG`) at the 3' cassette. Any additional site — introduced by a
    crRNA + scaffold junction that happens to spell CACCTGC or GCAGGTG across
    the boundary, or by a future code edit that breaks the wrapping — would
    cause PaqCI to cut in an unintended place and the assembly would fail.

    Called at the end of `build_fragment` and `Array.contiguous_insert` as a
    defensive guard. If this fires in the wild it means the user's specific
    crRNA is incompatible with its neighbors' sequences; the per-crRNA
    validator (which only sees the crRNA in isolation) can't catch these.

    Args:
        sequence: The full ordered DNA to scan.
        expected_fwd: Expected number of `CACCTGC` occurrences (default 1).
        expected_rev: Expected number of `GCAGGTG` occurrences (default 1).
        context: Human-readable label prepended to the error message
            (e.g. "position-3 fragment for gene Pax6").

    Raises:
        ValueError: If either count differs from the expected value. Message
            includes positions of every hit so the caller can pinpoint the
            reconstituted site (typically at a crRNA/scaffold or overhang/core
            boundary).
    """
    hits = scan_paqci_sites(sequence)
    fwd_hits = [i for i, o in hits if o == "fwd"]
    rev_hits = [i for i, o in hits if o == "rev"]
    if len(fwd_hits) == expected_fwd and len(rev_hits) == expected_rev:
        return
    prefix = f"{context}: " if context else ""
    raise ValueError(
        f"{prefix}unexpected PaqCI site count. "
        f"Expected {expected_fwd} CACCTGC (found {len(fwd_hits)} at {fwd_hits}) "
        f"and {expected_rev} GCAGGTG (found {len(rev_hits)} at {rev_hits}). "
        "This means a crRNA/scaffold/HDV boundary or barcode/stem-I boundary "
        "reconstitutes a PaqCI recognition site — the fragment will not "
        "assemble correctly. Pick a different crRNA or barcode."
    )


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
