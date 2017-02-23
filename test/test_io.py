from hw2skeleton import io
import pytest
import os

@pytest.mark.parametrize("filename,name,sequence", [
    ("prot-0004.fa", "d1flp__ 1.1.1.1.2", "SLEAAQKSNVTSSWAKASAAWGTAGPEFFMALFDAHDDVFAKFSGLFSGAAKGTVKNTPEMAAQAQSFKGLVSNWVDNLDNAGALEGQCKTFAANHKARGISAGQLEAAFKVLSGFMKSYGGDEGAWTAVAGALMGEIEPDM"),
])

def test_sequences(filename, name, sequence):
    filepath = os.path.join("sequences", filename)

    seq = io.read_FASTA(filepath)

    assert seq.name == name
    assert seq.sequence == sequence
