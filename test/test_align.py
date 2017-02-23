from hw2skeleton import align
from hw2skeleton import io
import os
import numpy as np
import math


def test_simple_align():
    seq1 = "ACV"
    seq2 = "ACV"

    score_matrix, key = align.read_score_matrix("BLOSUM50")

    max_score = align.SWalign(seq1,seq2,10,1,score_matrix,key)

    assert(max_score == 23)

def test_traceback():

    seq1 = "ACV"
    seq2 = "ACV"

    score_matrix, key = align.read_score_matrix("BLOSUM50")
    max_score, a1, a2 = align.SWalign_return_static(seq1,seq2,10,1,score_matrix,key)

    assert(a1 == "ACV")
    assert(a2 == "ACV")

