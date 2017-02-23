# A utility class to represent a FASTA sequence

class Sequence:
    """
    A simple class for an FASTA sequence
    """

    def __init__(self, name):
        self.name = name
        self.sequence = ""

    # Overload the __repr__ operator to make printing simpler.
    def __repr__(self):
        return self.name
