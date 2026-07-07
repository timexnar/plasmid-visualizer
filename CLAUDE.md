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
- [x] Stage 1 — Generalize input: `plasmid_io.load_record()` picks GenBank/FASTA/plain-text by extension; both scripts take an optional file path CLI arg (default: the sample file) and degrade gracefully with no annotations
- [x] Stage 2 — Restriction site detection: `restriction_sites.find_unique_cutters()` finds enzymes (from a curated common-cloning list) that cut exactly once, via `Bio.Restriction`. Listed in `parse_plasmid.py`; drawn as tick marks in their own outer ring in `plot_plasmid.py`, with a cluster-and-center declutter step (`declutter_label_angles()`) so tightly packed sites (e.g. an MCS) stay legible
- [x] Stage 3 — Insert matching: `insert_matching.find_insert_matches()` searches both strands and treats the plasmid as circular (matches spanning the origin are found and clipped for display). Both scripts gain `--insert`/`--insert-file`; `plot_plasmid.py` highlights matches as a distinct black arrow via `build_insert_features()`, reusing the existing feature pipeline
- [x] Stage 4 — Interactive map: `plot_plasmid_interactive.py` draws the same features/restriction sites/insert via Plotly `Barpolar`, exported as a standalone HTML file (embeds plotly.js) with hover tooltips and zoom/pan. Overlapping features are kept on separate rings via `assign_levels()`. Hover tooltips replace the need for the static map's label-declutter trick. `plasmid_map.html` is gitignored (a few MB from embedded plotly.js)
- [x] Stage 5 — Polish: `plot_plasmid.py` gained a dynamic legend (`draw_legend()`) and now decluttering feature labels too (not just restriction sites), via the same `declutter_label_angles()` with a wider gap for longer text (`FEATURE_MIN_LABEL_ANGLE_GAP`); fixed a sort-key bug where labels were ordered by feature start instead of midpoint, which could misplace whole-plasmid features like "source". `plot_plasmid_interactive.py` gained `--export-image` for a static PNG/SVG snapshot via kaleido.

All 5 planned stages are done. The tool now: parses/plots any GenBank, FASTA, or plain-text sequence file; detects and displays unique restriction cutters; highlights a user-supplied insert (both strands, circular-aware); offers both a static (matplotlib) and interactive (Plotly, hover/zoom, HTML/PNG/SVG export) circular map with legends and decluttered labels.

Possible future ideas (not currently planned): a Streamlit UI wrapping these scripts; smarter cross-ring label collision avoidance; support for multiple enzyme cuts (not just unique cutters) as an opt-in.
