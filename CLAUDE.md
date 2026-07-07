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
