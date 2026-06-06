"""Hammerhead and HDV ribozyme sequence assembly.

The HH ribozyme cleaves at the 5' end of the sgRNA spacer to give a precise 5' end.
Its Stem-I region is a 6-nt segment that base-pairs antiparallel with the first 6 nt
of the spacer (the crRNA). Stem-I is therefore *per crRNA* — it must be computed as
the reverse complement of `crRNA[:6]`.

The HDV ribozyme is constant. See `constants.HDV_RIBOZYME`.
"""

from __future__ import annotations

from sgrna_array.constants import (
    HDV_RIBOZYME,
    HH_CATALYTIC_CORE,
    HH_STEM_I_LEN,
)

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence, preserving input case."""
    return seq.translate(_COMPLEMENT)[::-1]


def build_hh_stem(crrna: str) -> str:
    """Build the 6-nt Stem-I segment for an HH ribozyme paired with this crRNA.

    Stem-I = reverse complement of the first 6 nt of the crRNA spacer. When transcribed,
    Stem-I pairs antiparallel with crRNA[:6] to form the 6-bp Stem I helix that positions
    the catalytic core to cleave precisely 5' of the spacer.

    Args:
        crrna: 20-nt crRNA spacer sequence. Must be at least HH_STEM_I_LEN nt.

    Returns:
        6-nt Stem-I DNA sequence (lowercase to match scaffold convention in dossier).

    Example:
        >>> build_hh_stem("ggtgacgcgtgtgtacctgg")
        'gtcacc'
    """
    if len(crrna) < HH_STEM_I_LEN:
        raise ValueError(
            f"crRNA must be at least {HH_STEM_I_LEN} nt to derive Stem-I; got {len(crrna)}"
        )
    return reverse_complement(crrna[:HH_STEM_I_LEN]).lower()


def build_full_hh(crrna: str) -> str:
    """Build the full HH ribozyme DNA encoding for a given crRNA: stem-I + catalytic core."""
    return build_hh_stem(crrna) + HH_CATALYTIC_CORE


def build_hdv() -> str:
    """Return the constant HDV ribozyme sequence."""
    return HDV_RIBOZYME
