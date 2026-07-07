"""
Load a plasmid file (GenBank, FASTA, or plain-text sequence) and draw an
interactive circular plasmid map, saved as a standalone HTML file.

Unlike the static map (plot_plasmid.py), we don't need to cram every
label around the circle - hovering over a feature, restriction site, or
insert shows its details in a tooltip, and the whole map can be zoomed
and panned. This uses Plotly's "Barpolar" chart type, which draws bars
in polar (angle + radius) coordinates - a natural fit for arcs around a
circular plasmid.
"""

import argparse

import plotly.graph_objects as go

from insert_matching import build_insert_features, find_insert_matches, load_insert_sequence
from plasmid_io import load_record
from restriction_sites import find_unique_cutters

DEFAULT_FILE = "addgene-plasmid-50005-sequence-513490.gbk"
DEFAULT_OUTPUT_HTML = "plasmid_map.html"

# "source" spans the entire plasmid and just restates its length; primer
# binding sites are numerous and mostly useful for cloning/sequencing
# work. Neither is useful for getting an overview of the plasmid.
IGNORED_FEATURE_TYPES = ("source", "primer_bind")

FEATURE_COLORS_BY_TYPE = {
    "CDS": "#f9d0c4",
    "promoter": "#ffe28a",
    "protein_bind": "#c9c9ff",
    "rep_origin": "#8affc1",
    "misc_feature": "#d1d1d1",
    "insert": "#000000",
}
DEFAULT_FEATURE_COLOR = "#7245dc"

BACKBONE_COLOR = "#e5e5e5"
RESTRICTION_SITE_COLOR = "firebrick"

BACKBONE_OUTER_R = 10
LEVEL_HEIGHT = 14  # radial distance between stacked feature levels
LEVEL_THICKNESS = 12  # radial thickness of each feature arc (leaves a small gap)
RESTRICTION_GAP = 10  # radial gap between the outermost feature level and the restriction ring


def assign_levels(intervals):
    """Greedily assign a level (0, 1, 2, ...) to each (start, end)
    interval so that overlapping intervals never share a level - e.g. two
    genes that overlap get drawn as two separate rings instead of on top
    of each other."""
    level_ends = []
    levels = []
    for start, end in intervals:
        placed = False
        for level, level_end in enumerate(level_ends):
            if start >= level_end:
                level_ends[level] = end
                levels.append(level)
                placed = True
                break
        if not placed:
            level_ends.append(end)
            levels.append(len(level_ends) - 1)
    return levels


def drawable_features(record):
    """Return (feature, start, end) for features we'll draw, excluding
    ignored types and zero-length spans."""
    drawable = []
    for feature in record.features:
        if feature.type in IGNORED_FEATURE_TYPES:
            continue
        start = int(feature.location.start)
        end = int(feature.location.end)
        if end > start:
            drawable.append((feature, start, end))
    return drawable


def add_backbone(fig, length):
    fig.add_trace(go.Barpolar(
        r=[BACKBONE_OUTER_R], theta=[180], width=[360], base=[0],
        marker_color=BACKBONE_COLOR, marker_line_color="black", marker_line_width=1,
        hovertext=[f"{length} bp"], hoverinfo="text",
        name="backbone", showlegend=False,
    ))


def add_features(fig, record):
    """Draw every feature as a Barpolar arc, grouped into one trace per
    feature type (so the legend shows one entry per type). Returns the
    highest level used, so the caller knows how far out it's safe to
    start drawing the next ring."""
    length = len(record.seq)
    features = drawable_features(record)
    levels = assign_levels([(start, end) for _, start, end in features])

    grouped = {}
    for (feature, start, end), level in zip(features, levels):
        span = end - start
        theta_center = 360 * (start + span / 2) / length
        theta_width = max(360 * span / length, 0.5)
        label = feature.qualifiers.get("label", ["(unnamed)"])[0]
        strand = feature.location.strand
        strand_symbol = "+" if strand == 1 else "-" if strand == -1 else "?"
        hover = f"{label}<br>{feature.type}, {start + 1}-{end} ({strand_symbol})"

        group = grouped.setdefault(
            feature.type, {"theta": [], "width": [], "base": [], "hover": []}
        )
        group["theta"].append(theta_center)
        group["width"].append(theta_width)
        group["base"].append(BACKBONE_OUTER_R + level * LEVEL_HEIGHT)
        group["hover"].append(hover)

    for feature_type, group in grouped.items():
        color = FEATURE_COLORS_BY_TYPE.get(feature_type, DEFAULT_FEATURE_COLOR)
        fig.add_trace(go.Barpolar(
            r=[LEVEL_THICKNESS] * len(group["theta"]),
            theta=group["theta"],
            width=group["width"],
            base=group["base"],
            marker_color=color,
            marker_line_color="black",
            marker_line_width=1,
            name=feature_type,
            hovertext=group["hover"],
            hoverinfo="text",
        ))

    return max(levels, default=-1)


def add_restriction_sites(fig, record, radius):
    cut_sites = find_unique_cutters(record.seq)
    if not cut_sites:
        return
    length = len(record.seq)
    thetas = [360 * position / length for _, position in cut_sites]
    hover = [f"{name} cuts at position {position}" for name, position in cut_sites]

    fig.add_trace(go.Scatterpolar(
        r=[radius] * len(cut_sites),
        theta=thetas,
        mode="markers",
        marker=dict(symbol="line-ns", size=12, color=RESTRICTION_SITE_COLOR, line_width=2),
        name="restriction site",
        hovertext=hover,
        hoverinfo="text",
    ))


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "input_file",
        nargs="?",
        default=DEFAULT_FILE,
        help=f"GenBank, FASTA, or plain-text sequence file (default: {DEFAULT_FILE})",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=DEFAULT_OUTPUT_HTML,
        help=f"Output HTML path (default: {DEFAULT_OUTPUT_HTML})",
    )
    insert_group = parser.add_mutually_exclusive_group()
    insert_group.add_argument(
        "--insert",
        help="A DNA sequence to search for and highlight on the map (e.g. ATGC...)",
    )
    insert_group.add_argument(
        "--insert-file",
        help="A file (FASTA or plain text) containing the insert sequence",
    )
    parser.add_argument(
        "--export-image",
        help="Also save a static snapshot (PNG or SVG, inferred from the file extension), "
        "e.g. --export-image plasmid_map.png",
    )
    return parser.parse_args()


def build_interactive_figure(record):
    """Build the interactive Plotly circular map figure for a loaded
    record (with any insert features already merged in). Used by both
    this script's CLI and app.py."""
    fig = go.Figure()
    add_backbone(fig, len(record.seq))
    max_level = add_features(fig, record)
    restriction_radius = BACKBONE_OUTER_R + (max_level + 1) * LEVEL_HEIGHT + RESTRICTION_GAP
    add_restriction_sites(fig, record, restriction_radius)

    fig.update_layout(
        title=record.name,
        polar=dict(
            bgcolor="white",
            radialaxis=dict(visible=False, range=[0, restriction_radius + 8]),
            angularaxis=dict(visible=False, rotation=90, direction="clockwise"),
        ),
        showlegend=True,
        width=750,
        height=750,
    )
    return fig


def main():
    args = parse_args()
    record = load_record(args.input_file)

    insert_seq = load_insert_sequence(args.insert, args.insert_file)
    if insert_seq is not None:
        matches = find_insert_matches(record.seq, insert_seq)
        record.features.extend(build_insert_features(matches, len(record.seq)))
        if not matches:
            print("Warning: no match found for the given insert sequence.")

    fig = build_interactive_figure(record)

    fig.write_html(args.output, include_plotlyjs=True)
    print(f"Saved interactive circular plasmid map to {args.output}")

    if args.export_image:
        fig.write_image(args.export_image)
        print(f"Saved static snapshot to {args.export_image}")


if __name__ == "__main__":
    main()
