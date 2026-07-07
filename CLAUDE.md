# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

A small plasmid/vector map visualizer built on Biopython. Given a GenBank, FASTA, or plain-text DNA sequence file, it parses features, finds unique restriction enzyme cut sites, can highlight a user-supplied insert sequence, and draws a circular plasmid map (static PNG or interactive HTML) — plus a Streamlit app (`app.py`) that ties all of this into one browser UI. `addgene-plasmid-50005-sequence-513490.gbk` (a pUC cloning vector from Addgene) is the bundled sample/default input.

## Setup and commands

The project uses a local virtualenv at `.venv` (already present, not committed).

```bash
# Create the venv (if not already present)
python3 -m venv .venv

# Install dependencies
.venv/bin/pip install -r requirements.txt

# Text summary: name, length, features, restriction sites, insert matches
.venv/bin/python parse_plasmid.py [file] [--insert SEQ | --insert-file FILE]

# Static circular map (PNG)
.venv/bin/python plot_plasmid.py [file] [-o out.png] [--insert SEQ | --insert-file FILE]

# Interactive circular map (standalone HTML, hover/zoom)
.venv/bin/python plot_plasmid_interactive.py [file] [-o out.html] [--export-image out.png]

# Browser app wrapping all of the above (upload a file, type an insert, pick a map style)
.venv/bin/streamlit run app.py
```

All three CLI scripts default to the bundled sample file if no path is given. Dependencies are pinned in `requirements.txt`.

## Architecture notes

- `parse_plasmid.py` reads a single hardcoded GenBank file path (`GENBANK_FILE` constant) via `Bio.SeqIO.read(..., "genbank")`, which assumes exactly one record in the file.
- Feature labels are pulled from the `label` qualifier (`feature.qualifiers.get("label", ...)`), since SnapGene/Addgene GenBank exports store the human-readable feature name there rather than in the standard `gene`/`product` qualifiers.
- `Bio.SeqFeature` locations are 0-based, half-open; the script converts the start to 1-based for display (`start = int(feature.location.start) + 1`) while leaving `end` as-is, matching GenBank's own 1-based inclusive convention.
- Strand is stored as `1`, `-1`, or `None` on `feature.location.strand`; `strand_symbol()` maps these to `+`, `-`, `?`.
- `plot_plasmid.py` renders a circular plasmid map using `dna_features_viewer`. A `PlasmidTranslator` subclass of `BiopythonTranslator` hides `primer_bind` features and colors features by type. Feature labels are drawn manually (not via the library's default layout) so each label radiates outward at its own feature's angle around the circle ("clock style") instead of being stacked in a row above the circle — see `draw_clock_style_labels()`. `build_static_figure(record)` builds and returns the matplotlib `Figure` without saving it, so both the CLI (`main()`) and `app.py` can reuse it.
- `plot_plasmid_interactive.py` renders the same map with Plotly's `Barpolar` chart type instead, as a standalone HTML file with hover tooltips (embeds plotly.js, so the output is a few MB and gitignored). `build_interactive_figure(record)` builds and returns the Plotly `Figure`, reused by both the CLI and `app.py`.
- `app.py` is a Streamlit app (`streamlit run app.py`) that wraps `load_record()`, `find_unique_cutters()`, `find_insert_matches()`/`build_insert_features()`, and both `build_*_figure()` functions behind a file uploader, an insert-sequence text box, and a static/interactive map toggle. An uploaded file is written to a temp file (via `tempfile.NamedTemporaryFile`) since `load_record()` expects a path, not a file object.

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

All 5 planned stages are done. The tool now: parses/plots any GenBank, FASTA, or plain-text sequence file; detects and displays unique restriction cutters; highlights a user-supplied insert (both strands, circular-aware); offers both a static (matplotlib) and interactive (Plotly, hover/zoom, HTML/PNG/SVG export) circular map with legends and decluttered labels; and wraps all of it in a Streamlit browser app (`app.py`, `streamlit run app.py`).

Possible future ideas (not currently planned): deploying `app.py` to Streamlit Community Cloud for a public shareable link; smarter cross-ring label collision avoidance; support for multiple enzyme cuts (not just unique cutters) as an opt-in.
