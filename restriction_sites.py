"""
Find restriction enzyme cut sites in a plasmid sequence using Bio.Restriction.

Only "unique cutters" (enzymes that cut the plasmid exactly once) are
returned by default, since those are the enzymes actually useful for
linearizing a plasmid or cutting out an insert for cloning - an enzyme
that cuts several times isn't usable that way.
"""

from Bio.Restriction import Analysis, RestrictionBatch

# A curated set of enzymes commonly found in multiple cloning sites,
# rather than Bio.Restriction's full ~600-enzyme catalog, which would
# bury the useful cutters under a lot of noise.
COMMON_CLONING_ENZYMES = [
    "EcoRI", "BamHI", "HindIII", "XhoI", "XbaI", "NotI", "SalI", "SmaI",
    "KpnI", "SacI", "PstI", "SphI", "NcoI", "NdeI", "SpeI", "ApaI",
    "NheI", "BglII", "ClaI", "EcoRV", "PvuII", "SacII", "AflII",
]


def find_unique_cutters(sequence, enzyme_names=COMMON_CLONING_ENZYMES):
    """Return a list of (enzyme_name, cut_position) for enzymes from
    `enzyme_names` that cut the given sequence exactly once.

    The sequence is treated as circular (as a plasmid is), and
    `cut_position` is 1-based, matching GenBank's convention.
    """
    batch = RestrictionBatch(enzyme_names)
    analysis = Analysis(batch, sequence, linear=False)
    results = analysis.full()

    unique_cuts = [
        (str(enzyme), positions[0])
        for enzyme, positions in results.items()
        if len(positions) == 1
    ]
    unique_cuts.sort(key=lambda item: item[1])
    return unique_cuts
