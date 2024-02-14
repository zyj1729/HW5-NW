# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")
    NW = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)
    assert (NW._ali_scores == np.array([[  0., -11., -12., -13.],
       [-11.,   5.,  -6.,  -7.],
       [-12.,  -6.,   4.,  -7.],
       [-13.,  -7.,  -1.,   5.],
       [-14.,  -8.,  -6.,   4.]])).all(), "The alignment matrix is incorrect"
    assert (NW._bt == np.array([['*', '<', '<', '<'],
       ['^', '`', '<', '<'],
       ['^', '^', '`', '<'],
       ['^', '^', '`', '`'],
       ['^', '^', '`', '`']], dtype='<U1')).all(), "The backtrace matrix is incorrect"
    assert (NW._gap == array([['*', 'g', 'g', 'g'],
       ['g', '*', 'g', 'g'],
       ['g', 'g', '*', 'g'],
       ['g', 'g', '*', '*'],
       ['g', 'g', '*', '*']], dtype='<U1')).all(), "The gap matrix is incorrect"
    pass
    
    

def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")
    NW = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)
    score, seqA_align, seqB_align = NW.align(seq3, seq4)
    assert score == 17, "Alignment score is incorrect for the test case"
    assert seqA_align == 'MAVHQLIRRP', "Aligned seqA is incorrect"
    assert seqB_align == 'M---QLIRHP', "Aligned seqB is incorrect"
    pass
    
    




