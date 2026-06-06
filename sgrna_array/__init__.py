"""sgrna_array — design ribozyme-flanked sgRNA DNA fragments for PaqCI Golden Gate assembly.

Pure-Python core for the sgRNA array webtool. Importable as a library, runnable via
the `sgrna-array` CLI entry point, and used by the Flask webapp under `webapp/`.

Public API (intended; implementations land in Phase 1):
    build_fragment(crRNA, position, array_size, **opts) -> Fragment
    build_array(fragments) -> Array
    validate_crrna(seq) -> list[ValidationResult]
"""

__version__ = "0.1.0"

from sgrna_array import constants  # noqa: F401  (re-export module for `sgrna_array.constants`)

__all__ = ["__version__", "constants"]
