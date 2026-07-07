# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

Small Python utility for parsing GenBank plasmid files (`.gb`/`.gbk`) using Biopython. The repo currently contains a single script, `parse_plasmid.py`, which reads `addgene-plasmid-50005-sequence-513490.gbk` (a pUC cloning vector exported from SnapGene/Addgene) and prints the plasmid name, length, and a table of annotated features (type, label, start/end position, strand).

## Setup and commands

The project uses a local virtualenv at `.venv` (already present, not committed).

```bash
# Create the venv (if not already present)
python3 -m venv .venv

# Install dependencies
.venv/bin/pip install -r requirements.txt

# Run the script
.venv/bin/python parse_plasmid.py
```

Dependencies are pinned in `requirements.txt` (`biopython`, `numpy`).

## Architecture notes

- `parse_plasmid.py` reads a single hardcoded GenBank file path (`GENBANK_FILE` constant) via `Bio.SeqIO.read(..., "genbank")`, which assumes exactly one record in the file.
- Feature labels are pulled from the `label` qualifier (`feature.qualifiers.get("label", ...)`), since SnapGene/Addgene GenBank exports store the human-readable feature name there rather than in the standard `gene`/`product` qualifiers.
- `Bio.SeqFeature` locations are 0-based, half-open; the script converts the start to 1-based for display (`start = int(feature.location.start) + 1`) while leaving `end` as-is, matching GenBank's own 1-based inclusive convention.
- Strand is stored as `1`, `-1`, or `None` on `feature.location.strand`; `strand_symbol()` maps these to `+`, `-`, `?`.
- `plot_plasmid.py` also reads the same hardcoded GenBank file and renders a circular plasmid map (`plasmid_map.png`) using `dna_features_viewer`. A `PlasmidTranslator` subclass of `BiopythonTranslator` hides `primer_bind` features and colors features by type. Feature labels are drawn manually (not via the library's default layout) so each label radiates outward at its own feature's angle around the circle ("clock style") instead of being stacked in a row above the circle — see `draw_clock_style_labels()`.

## Roadmap toward the long-term goal

Long-term goal: a general tool that takes a DNA sequence or GenBank file and draws an interactive circular plasmid map (origin, resistance markers, promoters, restriction sites, inserted constructs) — a simplified SnapGene, with shareable output.

Status:
- [x] Stage 0 — Parse a GenBank file and print features (`parse_plasmid.py`)
- [x] Stage 0 — Static circular map with colored, clock-style-labeled features (`plot_plasmid.py`)
- [ ] Stage 1 — Generalize input: accept any GenBank file via CLI arg, and accept a raw FASTA/plain DNA sequence (no annotations) as a fallback input
- [ ] Stage 2 — Restriction site detection via `Bio.Restriction`, drawn as tick marks/labels on the map
- [ ] Stage 3 — Let the user supply an "insert" sequence and highlight where it matches in the plasmid
- [ ] Stage 4 — Make the map interactive (leading option: redraw with Plotly for hover/zoom/pan and shareable HTML export; alternative: wrap in a Streamlit app)
- [ ] Stage 5 — Polish: legend, label-collision handling, PNG/SVG/HTML export options

Next up: Stage 1.
