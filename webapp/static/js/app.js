/* sgRNA Array Designer — minimal client-side enhancements:
   - Re-render the positional rows when the array size changes
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
      // Wire validation on any newly added inputs.
      rowsEl
        .querySelectorAll('input[data-validate-crrna]:not([data-wired])')
        .forEach(function (input) {
          attachCrrnaListeners(input);
        });
    });
  }

  function makeRowFromTemplate(template, pos) {
    const clone = template.cloneNode(true);
    clone.classList.remove("row-error");
    // Update labels and input names/IDs for the new position.
    clone.querySelector(".row-pos").textContent = "pos " + pos;
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
    return clone;
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
