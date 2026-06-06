"""Command-line entry point: `sgrna-array build inputs.csv -o out/`.

Reads a CSV of (gene, crRNA, position) rows, validates each crRNA, generates the array,
and writes FASTA + CSV order sheet + GenBank to the output directory.

Input CSV columns (header required):
    gene_label, crrna, position
Optional columns:
    array_size  (one of 4/6/8/10/12; defaults to max(position) snapped up to supported size)
    barcode     (12 nt; if omitted, auto-generated)
    use_cs1     (true/false; defaults to false)
"""

from __future__ import annotations

import argparse
import csv
import random
import sys
from pathlib import Path

from sgrna_array.assembler import build_array
from sgrna_array.barcode import generate_barcode
from sgrna_array.constants import SUPPORTED_ARRAY_SIZES
from sgrna_array.exporters import write_all
from sgrna_array.validator import has_errors, validate_crrna


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    return _cmd_build(args)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="sgrna-array",
        description=(
            "Design ribozyme-flanked sgRNA arrays for PaqCI Golden Gate assembly. "
            "Reads a CSV of crRNAs and writes FASTA / order-sheet / GenBank exports."
        ),
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    build = subparsers.add_parser("build", help="Build an array from a CSV of crRNAs.")
    build.add_argument("input_csv", type=Path, help="Path to input CSV.")
    build.add_argument(
        "-o", "--output-dir", type=Path, required=True, help="Directory for output files."
    )
    build.add_argument(
        "--basename", default="array", help="Output file basename (default: 'array')."
    )
    build.add_argument(
        "--seed",
        type=int,
        default=None,
        help="RNG seed for reproducible barcode generation.",
    )
    return parser


def _cmd_build(args: argparse.Namespace) -> int:
    if not args.input_csv.exists():
        print(f"error: input CSV not found: {args.input_csv}", file=sys.stderr)
        return 2

    rows = list(_read_rows(args.input_csv))
    if not rows:
        print(f"error: no data rows in {args.input_csv}", file=sys.stderr)
        return 2

    # Sort by position to ensure index order matches array order.
    rows.sort(key=lambda r: r["position"])

    # Validate all crRNAs up front; abort on any hard errors.
    fatal = False
    for row in rows:
        findings = validate_crrna(row["crrna"])
        for f in findings:
            print(
                f"{args.input_csv.name}:{row['_lineno']}: {row['crrna']!r}: {f}",
                file=sys.stderr,
            )
        if has_errors(findings):
            fatal = True
    if fatal:
        print("error: at least one crRNA failed validation; aborting", file=sys.stderr)
        return 1

    array_size = _resolve_array_size(rows)
    if array_size not in SUPPORTED_ARRAY_SIZES:
        print(
            f"error: inferred array_size={array_size} not in {SUPPORTED_ARRAY_SIZES}",
            file=sys.stderr,
        )
        return 2

    # Pad or truncate to the array size. (We require exact in the CSV; partial arrays will
    # error here.)
    if len(rows) != array_size:
        print(
            f"error: array_size={array_size} expects {array_size} crRNAs, got {len(rows)}",
            file=sys.stderr,
        )
        return 2

    crrnas = [r["crrna"] for r in rows]
    labels = [r["gene_label"] for r in rows]

    rng = random.Random(args.seed)
    barcode = rows[0].get("barcode") or generate_barcode(rng)

    array = build_array(
        crrnas=crrnas,
        gene_labels=labels,
        array_size=array_size,
        barcode=barcode,
        use_cs1=any(r.get("use_cs1", False) for r in rows),
    )

    paths = write_all(array, output_dir=args.output_dir, basename=args.basename)
    print(f"wrote {len(paths)} files to {args.output_dir}/:")
    for label, path in paths.items():
        print(f"  {label}: {path.name}")
    return 0


def _read_rows(path: Path):
    with path.open() as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            return
        required = {"gene_label", "crrna", "position"}
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(f"missing required CSV columns: {sorted(missing)}")
        for i, row in enumerate(reader, start=2):  # start=2: header is line 1
            yield {
                "_lineno": i,
                "gene_label": row.get("gene_label", "").strip(),
                "crrna": row["crrna"].strip(),
                "position": int(row["position"]),
                "array_size": (
                    int(row["array_size"]) if row.get("array_size") else None
                ),
                "barcode": row.get("barcode", "").strip() or None,
                "use_cs1": row.get("use_cs1", "").strip().lower() in {"true", "1", "yes"},
            }


def _resolve_array_size(rows: list[dict]) -> int:
    # Explicit array_size column wins if given consistently.
    explicit = {r["array_size"] for r in rows if r["array_size"] is not None}
    if len(explicit) == 1:
        return explicit.pop()
    if len(explicit) > 1:
        raise ValueError(f"conflicting array_size values in CSV: {sorted(explicit)}")
    # Otherwise pick the smallest supported size that fits.
    max_position = max(r["position"] for r in rows)
    for size in SUPPORTED_ARRAY_SIZES:
        if size >= max_position:
            return size
    raise ValueError(f"max position {max_position} exceeds largest supported array size")


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
