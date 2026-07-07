"""
Load a GenBank (.gb/.gbk) plasmid file and draw a circular plasmid map,
saved as a PNG image.

Feature labels are placed "clock style": each label sits out at the same
angle as its feature (like an hour hand pointing at that feature), instead
of all being stacked in a row above the circle.
"""

import math

import matplotlib.pyplot as plt
from Bio import SeqIO
from dna_features_viewer import BiopythonTranslator, CircularGraphicRecord
from dna_features_viewer.compute_features_levels import compute_features_levels

GENBANK_FILE = "addgene-plasmid-50005-sequence-513490.gbk"
OUTPUT_IMAGE = "plasmid_map.png"

# Primer binding sites are numerous and mostly useful for cloning/sequencing
# work, not for getting an overview of the plasmid, so we hide them here.
IGNORED_FEATURE_TYPES = ("primer_bind",)

FEATURE_COLORS_BY_TYPE = {
    "CDS": "#f9d0c4",
    "promoter": "#ffe28a",
    "protein_bind": "#c9c9ff",
    "rep_origin": "#8affc1",
    "misc_feature": "#d1d1d1",
}
DEFAULT_FEATURE_COLOR = "#7245dc"

LABEL_GAP = 0.3  # extra radius (in plot units) between the outermost arc and the labels


class PlasmidTranslator(BiopythonTranslator):
    ignored_features_types = IGNORED_FEATURE_TYPES

    def compute_feature_color(self, feature):
        return FEATURE_COLORS_BY_TYPE.get(feature.type, DEFAULT_FEATURE_COLOR)


def draw_clock_style_labels(ax, graphic_record, levels):
    """Draw each feature's label out at its own angle, like a clock hand,
    instead of using dna_features_viewer's default stacked-labels-on-top
    layout."""
    max_level = max(levels.values(), default=0)
    label_radius = (
        graphic_record.radius
        + (max_level + 1) * graphic_record.feature_level_height
        + LABEL_GAP
    )

    for feature, level in levels.items():
        mid_position = (feature.start + feature.end) / 2
        angle_deg = graphic_record.position_to_angle(mid_position)
        angle_rad = math.radians(angle_deg)
        cos_a, sin_a = math.cos(angle_rad), math.sin(angle_rad)

        # Where the feature's own arc ends, radius-wise.
        arc_radius = graphic_record.radius + (level + 0.5) * graphic_record.feature_level_height
        x0 = arc_radius * cos_a
        y0 = arc_radius * sin_a - graphic_record.radius

        # Where the label sits, further out along the same angle.
        x1 = label_radius * cos_a
        y1 = label_radius * sin_a - graphic_record.radius

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


def main():
    record = SeqIO.read(GENBANK_FILE, "genbank")

    # PlasmidTranslator turns Biopython SeqFeature objects into the
    # GraphicFeature objects dna_features_viewer needs to draw them,
    # applying our custom filtering and coloring rules.
    graphic_record = PlasmidTranslator().translate_record(
        record, record_class=CircularGraphicRecord
    )

    fig, ax = plt.subplots(1, figsize=(9, 9))
    graphic_record.initialize_ax(ax, draw_line=True, with_ruler=False)

    levels = compute_features_levels(graphic_record.features)
    for feature, level in levels.items():
        graphic_record.plot_feature(ax, feature, level)

    label_radius = draw_clock_style_labels(ax, graphic_record, levels)

    half_span = label_radius + 0.6
    ax.set_xlim(-half_span, half_span)
    ax.set_ylim(-graphic_record.radius - half_span, -graphic_record.radius + half_span)

    ax.set_title(record.name)
    fig.savefig(OUTPUT_IMAGE, bbox_inches="tight")
    print(f"Saved circular plasmid map to {OUTPUT_IMAGE}")


if __name__ == "__main__":
    main()
