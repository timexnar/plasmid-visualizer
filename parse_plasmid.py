"""
Load a GenBank (.gb/.gbk) plasmid file and print its basic info:
- plasmid name and length
- a list of all annotated features (type, name, start/end, strand)
"""

from Bio import SeqIO

GENBANK_FILE = "addgene-plasmid-50005-sequence-513490.gbk"


def strand_symbol(strand):
    if strand == 1:
        return "+"
    elif strand == -1:
        return "-"
    else:
        return "?"


def main():
    record = SeqIO.read(GENBANK_FILE, "genbank")

    print(f"Plasmid name: {record.name}")
    print(f"Length: {len(record.seq)} bp")
    print()
    print("Features:")

    for feature in record.features:
        feature_type = feature.type
        # Most SnapGene/Addgene files store the display name under "label"
        label = feature.qualifiers.get("label", ["(unnamed)"])[0]
        start = int(feature.location.start) + 1  # convert to 1-based for display
        end = int(feature.location.end)
        strand = strand_symbol(feature.location.strand)

        print(f"  - {feature_type:<15} {label:<25} {start}-{end} ({strand})")


if __name__ == "__main__":
    main()
