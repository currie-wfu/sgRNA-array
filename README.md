---
title: sgRNA Array Designer
emoji: 🧬
colorFrom: blue
colorTo: green
sdk: docker
app_port: 7860
pinned: false
license: mit
short_description: Multi-sgRNA arrays — ribozyme + PaqCI Golden Gate assembly
---

# sgRNA Array Designer

Design ribozyme-flanked sgRNA DNA fragments for one-pot PaqCI Golden Gate assembly
into multi-sgRNA arrays expressed from a Pol II intron.

Each generated fragment carries a [5' PaqCI cassette](https://www.neb.com/en-us/products/r0745-paqci) +
position-specific 4-nt overhang, an [HH ribozyme](https://en.wikipedia.org/wiki/Hammerhead_ribozyme)
with per-crRNA Stem-I, the crRNA spacer, the [DeWeirdt 2022](https://doi.org/10.1038/s41467-022-33024-2)
sgRNA scaffold, an [HDV ribozyme](https://en.wikipedia.org/wiki/Hepatitis_delta_virus_ribozyme),
and the matching 3' PaqCI cassette — ready to order from any oligo / gBlock vendor.

## Live tool

- **Web UI:** https://huggingface.co/spaces/currie-wfu/sgRNA_array
- **Landing / docs:** https://currie-wfu.github.io/sgRNA-array/

## Use locally

```bash
git clone https://github.com/currie-wfu/sgRNA-array.git
cd sgRNA-array

# Install
pip install -e .

# Run the test suite
python -m pytest tests/ -v

# Use the CLI for batch designs
python -m sgrna_array.cli build inputs.csv -o out/ --seed 42

# Run the Flask UI
pip install -e ".[web]"
flask --app webapp.app run --debug
```

## Validated design constraints

| Constraint | Value | Source |
|---|---|---|
| sgRNA scaffold | DeWeirdt 2022 (76 nt) | [10.1038/s41467-022-33024-2](https://doi.org/10.1038/s41467-022-33024-2) |
| HH catalytic core | 37 nt, constant | Host vector verified |
| HDV ribozyme | 68 nt, constant | Host vector verified |
| Type IIs enzyme | PaqCI (`CACCTGC`, 4-nt 5' overhang) | NEB R0745 |
| Overhang fidelity | 98% predicted (NEB Golden Gate Assembly Tool) | Potapov 2018 |
| Optional CS1 in Hairpin-1 loop | 22 nt | [10.1038/s41587-020-0470-y](https://doi.org/10.1038/s41587-020-0470-y) |
| Backbone overhangs (B5, B3) | `ACGG`, `GAGC` | Host vector PaqCI sites |

## Citing

A machine-readable citation lives in [`CITATION.cff`](CITATION.cff); GitHub renders it
in the right sidebar of the repo. A Zenodo-minted DOI will be added once a versioned
release is archived ([instructions](https://docs.github.com/en/repositories/archiving-a-github-repository/referencing-and-citing-content)).

```bibtex
@software{currie_sgrna_array_designer,
  author = {Currie, Joshua},
  title  = {sgRNA Array Designer},
  url    = {https://github.com/currie-wfu/sgRNA-array},
  year   = {2026},
  version = {0.1.0}
}
```

## License

MIT.
