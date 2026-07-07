"""
Load a plasmid file (GenBank, FASTA, or plain-text sequence) and draw a
circular plasmid map, saved as a PNG image.

Feature labels are placed "clock style": each label sits out at the same
angle as its feature (like an hour hand pointing at that feature), instead
of all being stacked in a row above the circle. Files with no annotations
(FASTA or plain text) simply produce an empty circle with no feature arcs.
"""

import argparse
import math

import matplotlib.pyplot as plt
from dna_features_viewer import BiopythonTranslator, CircularGraphicRecord
from dna_features_viewer.compute_features_levels import compute_features_levels
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from insert_matching import build_insert_features, find_insert_matches, load_insert_sequence
from plasmid_io import load_record
from restriction_sites import find_unique_cutters

DEFAULT_FILE = "addgene-plasmid-50005-sequence-513490.gbk"
DEFAULT_OUTPUT_IMAGE = "plasmid_map.png"

RESTRICTION_SITE_COLOR = "firebrick"
RESTRICTION_RING_GAP = 0.6  # radius gap between the feature labels and the restriction ring
RESTRICTION_TICK_LENGTH = 0.12

# Primer binding sites are numerous and mostly useful for cloning/sequencing
# work, not for getting an overview of the plasmid, so we hide them here.
IGNORED_FEATURE_TYPES = ("primer_bind",)

FEATURE_COLORS_BY_TYPE = {
    "CDS": "#f9d0c4",
    "promoter": "#ffe28a",
    "protein_bind": "#c9c9ff",
    "rep_origin": "#8affc1",
    "misc_feature": "#d1d1d1",
    "insert": "#000000",
}
DEFAULT_FEATURE_COLOR = "#7245dc"

LABEL_GAP = 0.3  # extra radius (in plot units) between the outermost arc and the labels


class PlasmidTranslator(BiopythonTranslator):
    ignored_features_types = IGNORED_FEATURE_TYPES

    def compute_feature_color(self, feature):
        return FEATURE_COLORS_BY_TYPE.get(feature.type, DEFAULT_FEATURE_COLOR)


# Minimum angular gap between neighboring labels. Feature names (e.g.
# "CAP binding site") are longer than enzyme names, so they need more
# angular room at the same label radius to avoid overlapping.
FEATURE_MIN_LABEL_ANGLE_GAP = 13
RESTRICTION_MIN_LABEL_ANGLE_GAP = 5


def declutter_label_angles(angles_deg, min_gap_deg):
    """Group tightly clustered angles and spread each group evenly around
    its own true center, so a busy region (like an MCS) fans out without
    dragging unrelated, well-separated labels along with it.

    `angles_deg` must be sorted highest to lowest (i.e. by ascending
    sequence position). Only the label position moves; the tick/arc stays
    at its true angle, connected to the label by a slightly bent leader
    line.
    """
    n = len(angles_deg)
    if n == 0:
        return []

    clusters = [[0]]
    for i in range(1, n):
        if angles_deg[i - 1] - angles_deg[i] < min_gap_deg:
            clusters[-1].append(i)
        else:
            clusters.append([i])

    adjusted = [None] * n
    for cluster in clusters:
        center = sum(angles_deg[i] for i in cluster) / len(cluster)
        k = len(cluster)
        for j, idx in enumerate(cluster):
            adjusted[idx] = center + min_gap_deg * ((k - 1) / 2 - j)
    return adjusted


def draw_clock_style_labels(ax, graphic_record, levels):
    """Draw each feature's label out at its own angle, like a clock hand,
    instead of using dna_features_viewer's default stacked-labels-on-top
    layout. Labels that end up too close together (e.g. inside an MCS)
    are spread apart via declutter_label_angles(), same as restriction
    site labels."""
    max_level = max(levels.values(), default=0)
    label_radius = (
        graphic_record.radius
        + (max_level + 1) * graphic_record.feature_level_height
        + LABEL_GAP
    )

    # Sort by the same midpoint used for the angle, so angles come out
    # highest-to-lowest, as declutter_label_angles() expects.
    ordered = sorted(levels.items(), key=lambda item: (item[0].start + item[0].end) / 2)
    tick_angles = [
        graphic_record.position_to_angle((feature.start + feature.end) / 2)
        for feature, _ in ordered
    ]
    label_angles = declutter_label_angles(tick_angles, FEATURE_MIN_LABEL_ANGLE_GAP)

    for (feature, level), tick_angle_deg, label_angle_deg in zip(
        ordered, tick_angles, label_angles
    ):
        tick_rad = math.radians(tick_angle_deg)
        # Where the feature's own arc ends, radius-wise.
        arc_radius = graphic_record.radius + (level + 0.5) * graphic_record.feature_level_height
        x0 = arc_radius * math.cos(tick_rad)
        y0 = arc_radius * math.sin(tick_rad) - graphic_record.radius

        # Where the label sits, further out, at its decluttered angle.
        label_rad = math.radians(label_angle_deg)
        x1 = label_radius * math.cos(label_rad)
        y1 = label_radius * math.sin(label_rad) - graphic_record.radius

        ax.plot([x0, x1], [y0, y1], color="0.5", linewidth=0.5, zorder=1, clip_on=False)

        # Labels on the right half of the circle read left-to-right away
        # from the circle; labels on the left half do the mirror image.
        ha = "left" if x1 >= 0 else "right"
        ax.text(
            x1,
            y1,
            feature.label,
            ha=ha,
            va="center",
            fontsize=8,
            clip_on=False,
        )

    return label_radius


def draw_restriction_sites(ax, graphic_record, cut_sites, inner_radius):
    """Draw a short tick mark and label for each unique restriction cut
    site, in its own ring beyond `inner_radius` so cut sites read as
    visually distinct from genetic features."""
    tick_outer_r = inner_radius + RESTRICTION_TICK_LENGTH
    label_r = tick_outer_r + 0.05

    # cut_sites is sorted by ascending position, which corresponds to
    # descending angle - exactly the order declutter_label_angles expects.
    tick_angles = [graphic_record.position_to_angle(position) for _, position in cut_sites]
    label_angles = declutter_label_angles(tick_angles, RESTRICTION_MIN_LABEL_ANGLE_GAP)

    for (enzyme_name, position), tick_angle_deg, label_angle_deg in zip(
        cut_sites, tick_angles, label_angles
    ):
        tick_rad = math.radians(tick_angle_deg)
        x0 = inner_radius * math.cos(tick_rad)
        y0 = inner_radius * math.sin(tick_rad) - graphic_record.radius

        label_rad = math.radians(label_angle_deg)
        x1 = tick_outer_r * math.cos(label_rad)
        y1 = tick_outer_r * math.sin(label_rad) - graphic_record.radius
        ax.plot(
            [x0, x1], [y0, y1],
            color=RESTRICTION_SITE_COLOR, linewidth=1, zorder=1, clip_on=False,
        )

        tx = label_r * math.cos(label_rad)
        ty = label_r * math.sin(label_rad) - graphic_record.radius
        ha = "left" if tx >= 0 else "right"
        ax.text(
            tx, ty,
            f"{enzyme_name} ({position})",
            ha=ha,
            va="center",
            fontsize=7,
            style="italic",
            color=RESTRICTION_SITE_COLOR,
            clip_on=False,
        )

    return label_r


def draw_legend(ax, record, cut_sites):
    """Add a legend listing only the feature types actually present on
    this map, plus a restriction site entry if any are drawn."""
    feature_types = sorted(
        {feature.type for feature in record.features if feature.type not in IGNORED_FEATURE_TYPES}
    )
    handles = [
        Patch(
            facecolor=FEATURE_COLORS_BY_TYPE.get(feature_type, DEFAULT_FEATURE_COLOR),
            edgecolor="black",
            linewidth=0.5,
            label=feature_type,
        )
        for feature_type in feature_types
    ]
    if cut_sites:
        handles.append(
            Line2D([0], [0], color=RESTRICTION_SITE_COLOR, label="restriction site")
        )
    ax.legend(handles=handles, loc="lower left", fontsize=8, frameon=False)


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
        default=DEFAULT_OUTPUT_IMAGE,
        help=f"Output PNG path (default: {DEFAULT_OUTPUT_IMAGE})",
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
    return parser.parse_args()


def build_static_figure(record):
    """Build the static matplotlib circular map figure for a loaded
    record (with any insert features already merged in). Used by both
    this script's CLI and app.py."""
    # PlasmidTranslator turns Biopython SeqFeature objects into the
    # GraphicFeature objects dna_features_viewer needs to draw them,
    # applying our custom filtering and coloring rules. Records with no
    # features (FASTA/plain text input) just produce an empty circle.
    graphic_record = PlasmidTranslator().translate_record(
        record, record_class=CircularGraphicRecord
    )

    fig, ax = plt.subplots(1, figsize=(9, 9))
    graphic_record.initialize_ax(ax, draw_line=True, with_ruler=False)

    levels = compute_features_levels(graphic_record.features)
    for feature, level in levels.items():
        graphic_record.plot_feature(ax, feature, level)

    label_radius = draw_clock_style_labels(ax, graphic_record, levels)

    cut_sites = find_unique_cutters(record.seq)
    restriction_radius = draw_restriction_sites(
        ax, graphic_record, cut_sites, inner_radius=label_radius + RESTRICTION_RING_GAP
    )

    draw_legend(ax, record, cut_sites)

    half_span = restriction_radius + 0.6
    ax.set_xlim(-half_span, half_span)
    ax.set_ylim(-graphic_record.radius - half_span, -graphic_record.radius + half_span)

    ax.set_title(record.name)
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

    fig = build_static_figure(record)
    fig.savefig(args.output, bbox_inches="tight")
    print(f"Saved circular plasmid map to {args.output}")


if __name__ == "__main__":
    main()
