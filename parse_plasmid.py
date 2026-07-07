"""
Load a plasmid file (GenBank, FASTA, or plain-text sequence) and print its
basic info:
- plasmid name and length
- a list of all annotated features (type, name, start/end, strand), if any
- a list of unique restriction enzyme cut sites
- where a user-supplied insert sequence matches the plasmid, if given
"""

import argparse

from insert_matching import find_insert_matches, load_insert_sequence
from plasmid_io import load_record
from restriction_sites import find_unique_cutters

DEFAULT_FILE = "addgene-plasmid-50005-sequence-513490.gbk"


def strand_symbol(strand):
    if strand == 1:
        return "+"
    elif strand == -1:
        return "-"
    else:
        return "?"


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "input_file",
        nargs="?",
        default=DEFAULT_FILE,
        help=f"GenBank, FASTA, or plain-text sequence file (default: {DEFAULT_FILE})",
    )
    insert_group = parser.add_mutually_exclusive_group()
    insert_group.add_argument(
        "--insert",
        help="A DNA sequence to search for in the plasmid (e.g. ATGC...)",
    )
    insert_group.add_argument(
        "--insert-file",
        help="A file (FASTA or plain text) containing the insert sequence",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    record = load_record(args.input_file)

    print(f"Plasmid name: {record.name}")
    print(f"Length: {len(record.seq)} bp")
    print()

    if record.features:
        print("Features:")
        for feature in record.features:
            feature_type = feature.type
            # Most SnapGene/Addgene files store the display name under "label"
            label = feature.qualifiers.get("label", ["(unnamed)"])[0]
            start = int(feature.location.start) + 1  # convert to 1-based for display
            end = int(feature.location.end)
            strand = strand_symbol(feature.location.strand)

            print(f"  - {feature_type:<15} {label:<25} {start}-{end} ({strand})")
    else:
        print("No annotated features found in this file.")

    print()
    unique_cutters = find_unique_cutters(record.seq)
    if unique_cutters:
        print("Unique restriction sites (cut exactly once):")
        for enzyme_name, position in unique_cutters:
            print(f"  - {enzyme_name:<10} cuts at position {position}")
    else:
        print("No unique restriction sites found among the common cloning enzymes checked.")

    insert_seq = load_insert_sequence(args.insert, args.insert_file)
    if insert_seq is not None:
        print()
        matches = find_insert_matches(record.seq, insert_seq)
        if matches:
            print("Insert matches:")
            for start, end, strand in matches:
                note = " (wraps around the origin)" if end > len(record.seq) else ""
                print(f"  - {start}-{end} ({strand_symbol(strand)}){note}")
        else:
            print("No match found for the given insert sequence.")


if __name__ == "__main__":
    main()
