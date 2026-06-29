"""Flask UI for the sgRNA array webtool.

Routes:
    GET  /                        — input form
    POST /design                  — validate inputs, build array, render results
    GET  /export/<design_id>/<fmt> — download FASTA / CSV / GenBank for a built array

The design state lives in a process-local dict keyed by a URL-safe random ID. This is
fine for development and small public deployments; for high-throughput / multi-worker
production, swap `_design_cache` for Redis or another shared store.

Run with:
    flask --app webapp.app run --debug
    # or
    python -m webapp.app
"""

from __future__ import annotations

import io
import random
import secrets
from typing import Any

from flask import Flask, abort, render_template, request, send_file

from sgrna_array.assembler import Array, build_array
from sgrna_array.barcode import generate_barcode, validate_barcode
from sgrna_array.constants import SUPPORTED_ARRAY_SIZES
from sgrna_array.exporters import to_csv_order_sheet, to_fasta, to_genbank
from sgrna_array.validator import Severity, validate_crrna


def create_app() -> Flask:
    """Application factory. Use this for testing; the module-level `app` is for `flask run`."""
    app = Flask(__name__, template_folder="templates", static_folder="static")
    app.config["MAX_CONTENT_LENGTH"] = 64 * 1024  # 64 KB easily covers a 12-row form

    # Process-local cache. Persists for the life of the Flask worker.
    design_cache: dict[str, Array] = {}
    app.config["DESIGN_CACHE"] = design_cache

    @app.route("/")
    def index() -> str:
        size = _coerce_size(request.args.get("size"), default=4)
        return render_template(
            "index.html",
            supported_sizes=SUPPORTED_ARRAY_SIZES,
            selected_size=size,
            form_errors=None,
            form_warnings=None,
            form_values=None,
        )

    @app.route("/design", methods=["POST"])
    def design() -> Any:
        array_size = _coerce_size(request.form.get("array_size"), default=None)
        if array_size is None:
            abort(400, "invalid array_size")

        labels, crrnas = _collect_rows(request.form, array_size)
        user_barcode = request.form.get("barcode", "").strip().upper() or None
        use_cs1 = request.form.get("use_cs1") == "on"

        errors, warnings = _validate_all(crrnas, user_barcode)
        if errors:
            return render_template(
                "index.html",
                supported_sizes=SUPPORTED_ARRAY_SIZES,
                selected_size=array_size,
                form_errors=errors,
                form_warnings=warnings,
                form_values={
                    "labels": labels,
                    "crrnas": crrnas,
                    "barcode": user_barcode or "",
                    "use_cs1": use_cs1,
                },
            )

        rng = random.Random()
        barcode = user_barcode or generate_barcode(rng)

        try:
            array = build_array(
                crrnas=crrnas,
                gene_labels=labels,
                array_size=array_size,
                barcode=barcode,
                use_cs1=use_cs1,
            )
        except ValueError as exc:
            abort(400, f"Array build failed: {exc}")

        design_id = secrets.token_urlsafe(16)
        design_cache[design_id] = array

        return render_template(
            "results.html",
            array=array,
            design_id=design_id,
            barcode=barcode,
            use_cs1=use_cs1,
            warnings=warnings,
            assembly_map_svg=_build_assembly_map_svg(array),
            contiguous_insert=array.contiguous_insert(),
        )

    @app.route("/export/<design_id>/<fmt>")
    def export(design_id: str, fmt: str) -> Any:
        array = design_cache.get(design_id)
        if array is None:
            abort(404, "Design not found or expired (restart the server clears the cache).")
        if fmt == "fasta":
            return _send_text(to_fasta(array), f"array_{design_id}.fasta", "text/plain")
        if fmt == "csv":
            return _send_text(to_csv_order_sheet(array), f"array_{design_id}.csv", "text/csv")
        if fmt == "genbank":
            return _send_text(to_genbank(array), f"array_{design_id}.gb", "text/plain")
        abort(404, "Unknown export format")

    return app


def _coerce_size(raw: str | None, default: int | None) -> int | None:
    try:
        value = int(raw or "")
    except ValueError:
        return default
    if value not in SUPPORTED_ARRAY_SIZES:
        return default
    return value


def _collect_rows(form: Any, array_size: int) -> tuple[list[str], list[str]]:
    labels: list[str] = []
    crrnas: list[str] = []
    for i in range(1, array_size + 1):
        labels.append(form.get(f"gene_label_{i}", "").strip())
        crrnas.append(form.get(f"crrna_{i}", "").strip().upper())
    return labels, crrnas


def _validate_all(
    crrnas: list[str], user_barcode: str | None
) -> tuple[dict[int, list[str]], dict[int, list[str]]]:
    """Return (errors, warnings) keyed by row index (1-based; 0 reserved for the barcode)."""
    errors: dict[int, list[str]] = {}
    warnings: dict[int, list[str]] = {}
    for i, crrna in enumerate(crrnas, start=1):
        results = validate_crrna(crrna)
        errs = [r.message for r in results if r.severity is Severity.ERROR]
        warns = [r.message for r in results if r.severity is Severity.WARNING]
        if errs:
            errors[i] = errs
        if warns:
            warnings[i] = warns
    if user_barcode:
        bc_results = validate_barcode(user_barcode)
        bc_errs = [r.message for r in bc_results if r.severity is Severity.ERROR]
        bc_warns = [r.message for r in bc_results if r.severity is Severity.WARNING]
        if bc_errs:
            errors[0] = bc_errs
        if bc_warns:
            warnings[0] = bc_warns
    return errors, warnings


def _send_text(content: str, filename: str, mimetype: str) -> Any:
    buf = io.BytesIO(content.encode("utf-8"))
    return send_file(buf, as_attachment=True, download_name=filename, mimetype=mimetype)


# ---------------------------------------------------------------------------
# Assembly map SVG
# ---------------------------------------------------------------------------


def _build_assembly_map_svg(array: Array) -> str:
    """Build a horizontal SVG of the array's fragment-and-junction layout.

    Each fragment is a rounded rectangle labeled with its gene + position; each junction
    is a narrow colored block whose hue is derived from the overhang sequence so the same
    overhang gets the same color across the diagram.
    """
    frag_w = 120
    junc_w = 48
    height = 80
    pad_x = 12
    pad_y = 16

    pieces: list[tuple[str, str, str]] = []
    # Each piece: (kind, label, color). Kind in {"junction", "fragment", "stitcher"}.
    pieces.append(("junction", array.fragments[0].left_overhang, _overhang_color(array.fragments[0].left_overhang)))
    for frag in array.fragments:
        pieces.append(("fragment", _fragment_label(frag), "#e7f0f7"))
        pieces.append(("junction", frag.right_overhang, _overhang_color(frag.right_overhang)))
    if array.stitcher is not None:
        pieces.append(("stitcher", f"stitcher\nstop-{array.stitcher.stop_size}", "#f7e7ef"))
        pieces.append(("junction", array.stitcher.right_overhang, _overhang_color(array.stitcher.right_overhang)))

    total_w = pad_x * 2 + sum(
        junc_w if kind == "junction" else frag_w for kind, _, _ in pieces
    )

    parts: list[str] = []
    parts.append(
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {total_w} {height + pad_y * 2}" '
        f'role="img" aria-label="Array assembly map">'
    )
    x = pad_x
    y = pad_y
    for kind, label, color in pieces:
        if kind == "junction":
            parts.append(
                f'<rect x="{x}" y="{y + 12}" width="{junc_w}" height="{height - 24}" '
                f'rx="4" ry="4" fill="{color}" stroke="#445" stroke-width="0.5"/>'
            )
            parts.append(
                f'<text x="{x + junc_w / 2}" y="{y + height / 2 + 4}" text-anchor="middle" '
                f'font-family="ui-monospace, monospace" font-size="12" fill="#222">{label}</text>'
            )
            x += junc_w
        else:
            parts.append(
                f'<rect x="{x}" y="{y}" width="{frag_w}" height="{height}" rx="8" ry="8" '
                f'fill="{color}" stroke="#244" stroke-width="1"/>'
            )
            lines = label.split("\n")
            for li, line in enumerate(lines):
                parts.append(
                    f'<text x="{x + frag_w / 2}" y="{y + height / 2 + (li - (len(lines) - 1) / 2) * 16 + 4}" '
                    f'text-anchor="middle" font-family="system-ui, sans-serif" font-size="13" '
                    f'fill="#114">{line}</text>'
                )
            x += frag_w
    parts.append("</svg>")
    return "".join(parts)


def _fragment_label(frag: Any) -> str:
    label = frag.gene_label or "(unnamed)"
    suffix = " · CS1" if frag.use_cs1 else ""
    return f"pos {frag.position}\n{label}{suffix}"


def _overhang_color(overhang: str) -> str:
    """Deterministic pastel color from a 4-nt overhang. Same overhang → same color."""
    h = (sum((ord(c) - 65) * (i + 1) for i, c in enumerate(overhang)) * 47) % 360
    return f"hsl({h}, 55%, 85%)"


# Module-level app for `flask --app webapp.app run` and `python -m webapp.app`.
app = create_app()


if __name__ == "__main__":  # pragma: no cover
    app.run(debug=True, port=5000)
