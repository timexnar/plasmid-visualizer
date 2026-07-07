"""
Find where a user-supplied "insert" sequence appears within a plasmid.

The plasmid is treated as circular (a match spanning the origin - the
start/end of the sequence - is still found), and both the forward strand
and its reverse complement are searched, since a cloned insert can end up
in either orientation.
"""

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature

INSERT_FEATURE_TYPE = "insert"


def find_insert_matches(plasmid_seq, insert_seq):
    """Return a list of (start, end, strand) 1-based matches of
    `insert_seq` within `plasmid_seq`. strand is 1 (forward) or -1
    (reverse complement). `end` may exceed the plasmid length for a
    match that wraps around the origin.
    """
    plasmid_seq = str(plasmid_seq).upper()
    insert_seq = str(insert_seq).upper()
    length = len(plasmid_seq)
    if not insert_seq:
        return []

    # Extend the sequence so a plain substring search can still find a
    # match that wraps past the last base back to the first.
    extended = plasmid_seq + plasmid_seq[: len(insert_seq) - 1]

    matches = []
    for strand, query in (
        (1, insert_seq),
        (-1, str(Seq(insert_seq).reverse_complement())),
    ):
        search_from = 0
        while True:
            index = extended.find(query, search_from)
            if index == -1 or index >= length:
                break
            start = index + 1
            end = index + len(query)
            matches.append((start, end, strand))
            search_from = index + 1

    return matches


def load_insert_sequence(raw_sequence=None, file_path=None):
    """Resolve the insert sequence from either a literal sequence string
    or a file path, loaded the same way plasmid input files are. Returns
    None if neither is given."""
    if raw_sequence:
        return raw_sequence
    if file_path:
        from plasmid_io import load_record

        return str(load_record(file_path).seq)
    return None


def build_insert_features(matches, plasmid_length):
    """Turn matches from find_insert_matches() into Biopython SeqFeatures,
    clipped to the plasmid length if a match wraps around the origin."""
    features = []
    for start, end, strand in matches:
        clipped_end = min(end, plasmid_length)
        location = FeatureLocation(start - 1, clipped_end, strand=strand)
        features.append(
            SeqFeature(location=location, type=INSERT_FEATURE_TYPE, qualifiers={"label": ["insert"]})
        )
    return features
