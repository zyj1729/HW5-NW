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
    test_NW = NeedlemanWunsch('substitution_matrices/BLOSUM62.mat', gap_open = -10, gap_extend = -1)
    alignment_score, seq1_alignment, seq2_alignment = test_NW.align(seq1, seq2)
    assert alignment_score == 4.0
    assert seq1_alignment == "MYQR"
    assert seq2_alignment == "M-QR"
    align_matrix = np.array([[  0., -np.inf, -np.inf, -np.inf],
 					[-np.inf,   5., -11., -13.],
 					[-np.inf, -12.,   4.,  -8.],
					[-np.inf, -12.,  -1.,   5.],
					[-np.inf, -14.,  -6.,   4.]])
    assert test_NW._align_matrix.all() == align_matrix.all()
    gapA_matrix = np.array([[-10., -np.inf, -np.inf, -np.inf],
							[-11., -12.,  -6.,  -7.],
 							[-12., -13., -14.,  -7.],
 							[-13., -14., -15., -12.],
 							[-14., -15., -16., -17.]])
    assert test_NW._gapA_matrix.all() == gapA_matrix.all()
    gapB_matrix = np.array([[-10., -11., -12., -13.],
 							[-np.inf, -12., -13., -14.],
 							[-np.inf,  -6., -14., -15.],
 							[-np.inf,  -7.,  -7., -16.],
 							[-np.inf,  -8.,  -8.,  -6.]])
    assert test_NW._gapB_matrix.all() == gapB_matrix.all()
    back = np.array([[-np.inf, -np.inf, -np.inf, -np.inf],
 					 [-np.inf,   0.,  -1.,  -1.],
 					 [-np.inf,   1.,   0.,  -1.],
					 [-np.inf,   1.,   0.,   0.],
					 [-np.inf,   1.,   0.,   0.]])
    assert test_NW._back.all() == back.all()
    

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
    test_NW = NeedlemanWunsch('substitution_matrices/BLOSUM62.mat', gap_open = -10, gap_extend = -1)
    alignment_score, seq3_alignment, seq4_alignment = test_NW.align(seq3, seq4)
    assert alignment_score == 17.0
    assert seq3_alignment == "MAVHQLIRRP"
    assert seq4_alignment == "M---QLIRHP"




