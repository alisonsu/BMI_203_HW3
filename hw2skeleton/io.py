import glob
import os
from Bio import SeqIO
from .utils import Sequence


def read_sequences(dir):
    """
    Read in all of the sequences from the given directory.

    Input: directory
    Output: list of ActiveSite instances
    """
    files = glob.glob(dir + '/*.fa')

    sequences = []
    # iterate over each .pdb file in the given directory
    for filepath in glob.iglob(os.path.join(dir, "*.fa")):

        sequences.append(read_FASTA(filepath))

    print("Read in %d sequences"%len(sequences))

    return sequences


def read_FASTA(filepath):
    """
    Read in a single active site given a PDB file

    Input: PDB file path
    Output: ActiveSite instance
    """
    basename = os.path.basename(filepath)
    name = os.path.splitext(basename)

    if name[1] != ".fa":
        raise IOError("%s is not a FASTA file"%filepath)

    r_num = 0

    # open pdb file
    for record in SeqIO.parse(filepath,"fasta"):
        name = record.description
        seq = Sequence(name)
        seq.sequence = str(record.seq)

    return(seq)
