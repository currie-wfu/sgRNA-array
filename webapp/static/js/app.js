/* sgRNA Array Designer — minimal client-side enhancements:
   - Re-render the positional rows when the array size changes
   - Drag-to-reorder rows via a grip handle (HTML5 native DnD)
   - Live crRNA validation hints (length, ACGT, PaqCI sites)
   - Copy-to-clipboard for fragment sequences
   No build step. Vanilla DOM only. */

(function () {
  "use strict";

  const PAQCI_FWD = "CACCTGC";
  const PAQCI_REV = "GCAGGTG";
  const CRRNA_LEN = 20;

  document.addEventListener("DOMContentLoaded", function () {
    wireArraySizeSelector();
    wireCrrnaInputs();
    wireDragAndDrop();
    wireCopyButtons();
    // crRNA inputs may have been server-rendered with values; validate them on load.
    document.querySelectorAll("[data-validate-crrna]").forEach(validateCrrna);
  });

  /* -------------------------------------------------------------------- */
  /* Array-size dropdown rebuilds the positional rows in-place.            */
  /* -------------------------------------------------------------------- */
  function wireArraySizeSelector() {
    const select = document.getElementById("array_size");
    const rowsEl = document.getElementById("rows");
    if (!select || !rowsEl) return;

    select.addEventListener("change", function () {
      const target = parseInt(select.value, 10);
      if (!Number.isFinite(target)) return;
      const existingRows = rowsEl.querySelectorAll(".row");
      const current = existingRows.length;

      if (target === current) return;

      if (target > current) {
        // Add rows for new positions, copying the structure of the first row.
        const template = existingRows[0];
        for (let pos = current + 1; pos <= target; pos++) {
          rowsEl.appendChild(makeRowFromTemplate(template, pos));
        }
      } else {
        // Trim extra rows.
        for (let i = current - 1; i >= target; i--) {
          rowsEl.removeChild(existingRows[i]);
        }
      }
      // Wire validation + DnD on any newly added rows.
      rowsEl
        .querySelectorAll('input[data-validate-crrna]:not([data-wired])')
        .forEach(attachCrrnaListeners);
      rowsEl
        .querySelectorAll('.row:not([data-dnd-wired])')
        .forEach(attachRowDragEvents);
      rowsEl
        .querySelectorAll('.row-drag-handle:not([data-dnd-wired])')
        .forEach(attachDragHandle);
      renumberRows();
    });
  }

  function makeRowFromTemplate(template, pos) {
    const clone = template.cloneNode(true);
    clone.classList.remove("row-error", "row-dragging", "row-drop-above", "row-drop-below");
    clone.removeAttribute("data-dnd-wired");
    // Update labels and input names/IDs for the new position.
    clone.querySelector(".row-pos").textContent = "pos " + pos;
    clone.querySelectorAll("input").forEach(function (input) {
      const name = input.getAttribute("name") || "";
      const newName = name.replace(/_\d+$/, "_" + pos);
      input.setAttribute("name", newName);
      input.value = "";
      input.classList.remove("row-input-warn");
      input.removeAttribute("data-wired");
    });
    // Remove any inherited error/warning lists.
    clone.querySelectorAll(".row-error-msgs, .row-warning-msgs").forEach(function (el) {
      el.remove();
    });
    // Clear DnD-wired marker on the handle so it gets re-wired.
    const handle = clone.querySelector(".row-drag-handle");
    if (handle) handle.removeAttribute("data-dnd-wired");
    return clone;
  }

  /* -------------------------------------------------------------------- */
  /* Drag-to-reorder rows. The handle toggles draggable on its parent     */
  /* row only while the handle is pressed; this keeps text-selection in   */
  /* the inputs working normally.                                          */
  /* -------------------------------------------------------------------- */
  function wireDragAndDrop() {
    document.querySelectorAll(".row").forEach(attachRowDragEvents);
    document.querySelectorAll(".row-drag-handle").forEach(attachDragHandle);
  }

  function attachDragHandle(handle) {
    if (handle.dataset.dndWired === "true") return;
    handle.dataset.dndWired = "true";
    const row = handle.closest(".row");
    if (!row) return;

    const enable = function () { row.setAttribute("draggable", "true"); };
    const disable = function () { row.removeAttribute("draggable"); };

    handle.addEventListener("mousedown", enable);
    handle.addEventListener("mouseup", disable);
    handle.addEventListener("mouseleave", disable);
    // Keyboard parity: allow Space/Enter on a focused handle to enable
    // drag, then standard DnD takes over once the user begins.
    handle.addEventListener("keydown", function (e) {
      if (e.key === " " || e.key === "Enter") {
        e.preventDefault();
        enable();
      }
    });
  }

  function attachRowDragEvents(row) {
    if (row.dataset.dndWired === "true") return;
    row.dataset.dndWired = "true";

    row.addEventListener("dragstart", function (e) {
      row.classList.add("row-dragging");
      e.dataTransfer.effectAllowed = "move";
      // Setting any data is required for Firefox to dispatch drag events.
      e.dataTransfer.setData("text/plain", "row");
    });

    row.addEventListener("dragend", function () {
      row.classList.remove("row-dragging");
      row.removeAttribute("draggable");
      document.querySelectorAll(".row-drop-above, .row-drop-below").forEach(function (el) {
        el.classList.remove("row-drop-above", "row-drop-below");
      });
      renumberRows();
    });

    row.addEventListener("dragover", function (e) {
      const dragging = document.querySelector(".row-dragging");
      if (!dragging || dragging === row) return;
      e.preventDefault();
      e.dataTransfer.dropEffect = "move";

      const rect = row.getBoundingClientRect();
      const midY = rect.top + rect.height / 2;
      const dropAbove = e.clientY < midY;

      row.classList.toggle("row-drop-above", dropAbove);
      row.classList.toggle("row-drop-below", !dropAbove);
    });

    row.addEventListener("dragleave", function () {
      row.classList.remove("row-drop-above", "row-drop-below");
    });

    row.addEventListener("drop", function (e) {
      e.preventDefault();
      const dragging = document.querySelector(".row-dragging");
      if (!dragging || dragging === row) return;

      const rect = row.getBoundingClientRect();
      const dropAbove = e.clientY < rect.top + rect.height / 2;
      const parent = row.parentNode;
      if (dropAbove) {
        parent.insertBefore(dragging, row);
      } else {
        parent.insertBefore(dragging, row.nextSibling);
      }
      row.classList.remove("row-drop-above", "row-drop-below");
    });
  }

  function renumberRows() {
    const rows = document.querySelectorAll("#rows .row");
    rows.forEach(function (row, idx) {
      const pos = idx + 1;
      const posLabel = row.querySelector(".row-pos");
      if (posLabel) posLabel.innerHTML = "pos&nbsp;" + pos;
      row.querySelectorAll("input").forEach(function (input) {
        const name = input.getAttribute("name") || "";
        const newName = name.replace(/_\d+$/, "_" + pos);
        if (newName !== name) input.setAttribute("name", newName);
      });
    });
  }

  /* -------------------------------------------------------------------- */
  /* crRNA input validation: length, alphabet, internal PaqCI sites.       */
  /* -------------------------------------------------------------------- */
  function wireCrrnaInputs() {
    document
      .querySelectorAll("input[data-validate-crrna]")
      .forEach(attachCrrnaListeners);
  }

  function attachCrrnaListeners(input) {
    if (input.dataset.wired === "true") return;
    input.dataset.wired = "true";
    input.addEventListener("input", function () {
      input.value = input.value.replace(/\s+/g, "").toUpperCase();
      validateCrrna(input);
    });
    input.addEventListener("blur", function () {
      validateCrrna(input);
    });
  }

  function validateCrrna(input) {
    const value = input.value.trim();
    input.classList.remove("row-input-warn");
    input.setCustomValidity("");
    if (value.length === 0) return;

    if (!/^[ACGT]+$/.test(value)) {
      input.setCustomValidity("crRNA must be ACGT only");
      return;
    }
    if (value.length !== CRRNA_LEN) {
      input.setCustomValidity("crRNA must be exactly " + CRRNA_LEN + " nt");
      return;
    }
    if (value.includes(PAQCI_FWD) || value.includes(PAQCI_REV)) {
      input.setCustomValidity("crRNA contains internal PaqCI recognition site");
      return;
    }
    // GC sanity (soft warning, not invalid).
    const gc = countGC(value) / value.length;
    if (gc < 0.25 || gc > 0.75) {
      input.classList.add("row-input-warn");
    }
  }

  function countGC(s) {
    let count = 0;
    for (let i = 0; i < s.length; i++) {
      if (s[i] === "G" || s[i] === "C") count++;
    }
    return count;
  }

  /* -------------------------------------------------------------------- */
  /* Copy-to-clipboard for fragment sequences on the results page.        */
  /* -------------------------------------------------------------------- */
  function wireCopyButtons() {
    document.querySelectorAll(".btn-copy").forEach(function (btn) {
      btn.addEventListener("click", function () {
        const targetId = btn.getAttribute("data-copy-target");
        const src = document.getElementById(targetId);
        if (!src) return;
        src.select();
        try {
          document.execCommand("copy");
          flashCopied(btn);
        } catch (e) {
          // Fallback to Clipboard API (newer browsers).
          if (navigator.clipboard) {
            navigator.clipboard.writeText(src.value).then(function () {
              flashCopied(btn);
            });
          }
        }
        window.getSelection().removeAllRanges();
      });
    });
  }

  function flashCopied(btn) {
    const originalText = btn.textContent;
    btn.classList.add("copied");
    btn.textContent = "Copied!";
    setTimeout(function () {
      btn.classList.remove("copied");
      btn.textContent = originalText;
    }, 1200);
  }
})();
