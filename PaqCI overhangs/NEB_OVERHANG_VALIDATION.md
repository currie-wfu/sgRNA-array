# NEB Golden Gate Assembly Tool — PaqCI overhang validation

This file documents the candidate 13-overhang set for the sgRNA array webtool and walks through validating it with NEB's free online tool. Once the validation is done, paste the result back in chat so the set can be locked into `sgrna_array/constants.py`.

---

## Background

A 12-position sgRNA array assembled by one-pot PaqCI Golden Gate needs **13 unique 4-nt overhangs**: one for the backbone-to-position-1 junction (`B5`), eleven inter-position junctions (`J1`…`J11`), and one for the position-12-to-backbone junction (`B3`).

**B5 and B3 are fixed by the existing host vector** (`1103p_SceITol2_14xUAS_CRISPR2-0ShhmNeonGreen-T2A-Cas9.gb`). PaqCI digestion of that vector creates two specific 4-nt overhangs that any new array MUST ligate into. The NEB tool must score the full 13-junction set including these locked values.

**Target:** ≥95% predicted fidelity for the full 13-junction set. Lower scores mean a meaningful fraction of clones will misassemble, dropping below the 25–30% correct-assembly screening threshold once colony screening losses compound.

---

## Candidate overhang set (B5/B3 locked, J1–J11 tunable)

| Junction | Role | Overhang (5'→3') | Status |
|----------|------|------------------|--------|
| **B5**  | backbone → position 1 | **`ACGG`** | 🔒 locked by vector (PaqCI rev cut at vector pos 1907-1910) |
| J1  | position 1 → position 2 | `AATG` | candidate |
| J2  | position 2 → position 3 | `TTCT` | candidate |
| J3  | position 3 → position 4 | `GGAG` | candidate |
| J4  | position 4 → position 5 (also = stop-at-4 stitcher overhang) | `AAGC` | candidate |
| J5  | position 5 → position 6 | `GCAA` | candidate |
| J6  | position 6 → position 7 (also = stop-at-6 stitcher overhang) | `CGAA` | candidate |
| J7  | position 7 → position 8 | `GTCT` | candidate |
| J8  | position 8 → position 9 (also = stop-at-8 stitcher overhang) | `ATCC` | candidate |
| J9  | position 9 → position 10 | `TGAC` | candidate |
| J10 | position 10 → position 11 (also = stop-at-10 stitcher overhang) | `CAGT` | candidate |
| J11 | position 11 → position 12 | `ACTC` | candidate |
| **B3**  | position 12 → backbone | **`GAGC`** | 🔒 locked by vector (PaqCI fwd cut at vector pos 2144-2147) |

### Why B5/B3 are locked
The host vector's two PaqCI cut sites release a fragment with these exact 4-nt overhangs. New array fragments must terminate with them to ligate into the backbone. The only way to change B5/B3 would be to re-clone the vector with different PaqCI cassette overhangs — out of scope for this tool.

### Design rationale for J1–J11
- All 11 are non-palindromic (no self-annealing).
- GC content per overhang: 25–75% (no `AAAA`/`TTTT` weakness, no `GGGG`/`CCCC` mispairing).
- Pairwise Hamming distance ≥2 within the set.
- None overlap the PaqCI recognition motif (`CACCTGC` or `GCAGGTG`) in flanking context.
- None coincide with `ACGG` or `GAGC` (the locked B5/B3).
- Drawn from published high-fidelity 4-nt overhang data (Potapov et al. 2018, ACS Synth Biol).

---

## How to run the NEB tool

1. **Open the tool:** go to <https://goldengate.neb.com>. No login required.
2. **Click "Manual Design"** (or the equivalent option for entering your own overhangs).
3. **Set the Type IIs enzyme to PaqCI.** This is in a dropdown near the top of the design form.
4. **Set the assembly conditions** (T4 DNA ligase, 37 °C, overnight is standard; if the tool offers presets, use the "high-fidelity, long incubation" preset).
5. **Enter the 13 overhangs in order:** `ACGG`, `AATG`, `TTCT`, `GGAG`, `AAGC`, `GCAA`, `CGAA`, `GTCT`, `ATCC`, `TGAC`, `CAGT`, `ACTC`, `GAGC`. The tool may ask for "fusion sites" or "junctions" — same thing.
6. **Click "Calculate Fidelity"** (or "Predict assembly" / "Evaluate").
7. **Record the result.** Looking for:
   - **Overall predicted fidelity (%)** — top-line number, ideally ≥95.
   - **Per-junction fidelity breakdown** if shown — flags any single weak junction.
   - **Misligation matrix** if shown — flags specific overhang pairs that cross-react.
   - **Tool-suggested improvements** — NEB sometimes offers J1–J11 swap candidates. (Note: don't accept swaps for B5 or B3 — those are locked by the vector.)

---

## What to paste back

Reply in chat with whatever the tool reports. Minimum useful info:

- Predicted overall fidelity, e.g. `97.4%`
- Any per-junction scores below ~95%
- Any J1–J11 swaps the tool recommended (NEB may not respect that B5/B3 are locked — ignore any suggestions to change those two)
- Optionally: a screenshot of the result page

If the tool reports the set is below 95%, prefer adjustments to J1–J11 over B5/B3 (which can't change without re-cloning the vector).

---

## If the tool can't be run

The same fidelity calculation can be run offline using Potapov et al.'s open ligation-bias dataset and `GGTools` at <https://github.com/potapovneb/golden-gate>. NEB's web tool is faster for a one-shot check; the offline path is only useful if we end up needing to programmatically search a larger overhang space.
