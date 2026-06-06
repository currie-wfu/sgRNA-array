"""Exporters for assembled arrays: FASTA, CSV order sheet, and annotated GenBank.

GenBank export uses BioPython's `SeqRecord` / `SeqFeature` and annotates the structural
modules of each fragment (PaqCI cassette, ribozymes, crRNA, scaffold, barcode) so the
output is reviewable in Snapgene / ApE / Benchling.
"""

from __future__ import annotations

import csv
import io
from pathlib import Path

from sgrna_array.assembler import Array, Fragment, Stitcher


def to_fasta(array: Array) -> str:
    """Return a FASTA-formatted string with one record per fragment + optional stitcher."""
    out = io.StringIO()
    for frag in array.fragments:
        header = _fasta_header_for_fragment(frag)
        out.write(f">{header}\n{frag.ordered_dna}\n")
    if array.stitcher is not None:
        header = _fasta_header_for_stitcher(array.stitcher)
        out.write(f">{header}\n{array.stitcher.ordered_dna}\n")
    return out.getvalue()


def to_csv_order_sheet(array: Array) -> str:
    """Return a CSV string suitable for upload to an oligo / gBlock vendor.

    Columns: piece_id, type, position, gene_label, length, sequence, notes
    """
    out = io.StringIO()
    writer = csv.writer(out)
    writer.writerow(
        ["piece_id", "type", "position", "gene_label", "length", "sequence", "notes"]
    )
    for frag in array.fragments:
        writer.writerow(
            [
                _piece_id(frag),
                "fragment",
                frag.position,
                frag.gene_label,
                len(frag.ordered_dna),
                frag.ordered_dna,
                _fragment_notes(frag, array.array_size),
            ]
        )
    if array.stitcher is not None:
        writer.writerow(
            [
                f"stitcher_stop_{array.stitcher.stop_size}",
                "stitcher",
                "",
                "",
                len(array.stitcher.ordered_dna),
                array.stitcher.ordered_dna,
                f"Terminates array at position {array.stitcher.stop_size}",
            ]
        )
    return out.getvalue()


def to_genbank(array: Array) -> str:
    """Return a multi-record GenBank string with annotated features per fragment.

    Each fragment becomes one GenBank record annotating: 5' pad, PaqCI fwd recog, left
    overhang, HH stem-I, HH core, crRNA, scaffold (with optional CS1 sub-feature), HDV,
    right overhang, PaqCI rev recog, 3' pad. The barcode is annotated on position-1.
    """
    # Lazy import so the rest of the library doesn't require BioPython at import time.
    from Bio.Seq import Seq
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord

    records: list[SeqRecord] = []
    for frag in array.fragments:
        record = SeqRecord(
            Seq(frag.ordered_dna.upper()),
            id=_piece_id(frag),
            name=_piece_id(frag)[:16],  # GenBank LOCUS limit
            description=(
                f"sgRNA array fragment position {frag.position}/{array.array_size}; "
                f"target {frag.gene_label or 'unspecified'}"
            ),
            annotations={"molecule_type": "DNA"},
        )
        record.features = list(_annotate_fragment(frag))
        records.append(record)

    if array.stitcher is not None:
        s = array.stitcher
        record = SeqRecord(
            Seq(s.ordered_dna.upper()),
            id=f"stitcher_stop_{s.stop_size}",
            name=f"stitch_stop{s.stop_size}"[:16],
            description=f"Stitcher oligo terminating array at position {s.stop_size}",
            annotations={"molecule_type": "DNA"},
        )
        records.append(record)

    out = io.StringIO()
    from Bio import SeqIO

    SeqIO.write(records, out, "genbank")
    return out.getvalue()


def write_all(array: Array, output_dir: str | Path, basename: str = "array") -> dict[str, Path]:
    """Write FASTA, CSV order sheet, and GenBank to `output_dir/{basename}.{ext}`."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "fasta": output_dir / f"{basename}.fasta",
        "csv": output_dir / f"{basename}.csv",
        "genbank": output_dir / f"{basename}.gb",
    }
    paths["fasta"].write_text(to_fasta(array))
    paths["csv"].write_text(to_csv_order_sheet(array))
    paths["genbank"].write_text(to_genbank(array))
    return paths


def _fasta_header_for_fragment(frag: Fragment) -> str:
    label = frag.gene_label.replace(" ", "_") or "unnamed"
    cs1 = "_cs1" if frag.use_cs1 else ""
    bc = f"_bc{frag.barcode}" if frag.barcode else ""
    return (
        f"pos{frag.position}_{label}{cs1}{bc} "
        f"L={frag.left_overhang} R={frag.right_overhang} "
        f"len={len(frag.ordered_dna)}"
    )


def _fasta_header_for_stitcher(stitcher: Stitcher) -> str:
    return (
        f"stitcher_stop{stitcher.stop_size} "
        f"L={stitcher.left_overhang} R={stitcher.right_overhang} "
        f"len={len(stitcher.ordered_dna)}"
    )


def _piece_id(frag: Fragment) -> str:
    label = frag.gene_label.replace(" ", "_") or "unnamed"
    return f"pos{frag.position}_{label}"


def _fragment_notes(frag: Fragment, array_size: int) -> str:
    notes = [f"pos {frag.position} of {array_size}"]
    if frag.barcode:
        notes.append(f"barcode={frag.barcode}")
    if frag.use_cs1:
        notes.append("CS1 in scaffold Hairpin-1 loop")
    return "; ".join(notes)


def _annotate_fragment(frag: Fragment):  # type: ignore[no-untyped-def]
    """Generate BioPython SeqFeature objects for the structural modules of a fragment.

    Implementation note: feature coordinate calculation is deferred to Task #4
    (golden-fragment tests) — for now this is a stub that returns no features so that
    GenBank export at least succeeds with bare sequences. A follow-up implementation
    will compute precise FeatureLocations for each module.
    """
    # TODO(Task #4 / sequence_design.md): annotate pad/recog/overhang/HH/crRNA/scaffold/HDV
    # using FeatureLocation with exact coordinates derived from constants module lengths.
    return iter(())
