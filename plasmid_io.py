"""
Shared helper for loading a DNA sequence from a file, used by both
parse_plasmid.py and plot_plasmid.py.

Three kinds of input are supported:
- GenBank (.gb, .gbk, .genbank): comes with annotated features
  (promoters, genes, origins, etc.)
- FASTA (.fasta, .fa, .fna): just a sequence and a header, no annotations
- Plain text: any other file is treated as raw sequence text (e.g. a
  sequence pasted into a .txt file), also with no annotations
"""

from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

GENBANK_EXTENSIONS = {".gb", ".gbk", ".genbank"}
FASTA_EXTENSIONS = {".fasta", ".fa", ".fna"}


def load_record(path):
    """Load a Biopython SeqRecord from a GenBank, FASTA, or plain-text file.

    Records loaded from FASTA or plain text have an empty ``features``
    list, since those formats carry no annotations.
    """
    suffix = Path(path).suffix.lower()

    if suffix in GENBANK_EXTENSIONS:
        return SeqIO.read(path, "genbank")

    if suffix in FASTA_EXTENSIONS:
        return SeqIO.read(path, "fasta")

    raw_text = Path(path).read_text()
    sequence = "".join(raw_text.split())
    name = Path(path).stem
    return SeqRecord(Seq(sequence), id=name, name=name)
