# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    # TODO Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    test_NW = NeedlemanWunsch('substitution_matrices/BLOSUM62.mat', gap_open = -10, gap_extend = -1)
    gg, _, _ = test_NW.align(hs_seq, gg_seq)
    mm, _, _ = test_NW.align(hs_seq, mm_seq)
    br, _, _ = test_NW.align(hs_seq, br_seq)
    tt, _, _ = test_NW.align(hs_seq, tt_seq)
    scores = {"Gallus gallus (red junglefowl)": gg, 
    	"Mus musculus (mouse)": mm, 
    	"Balaeniceps rex (shoebill)": br, 
    	"Tursiops truncatus (bottlenose dolphin)": tt}

    # TODO print all of the alignment score between each species BRD2 and human BRD2
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    print("ALIGNMENT SCORES TO HUMAN BRD2")
    scores = sorted(scores.items(), key = lambda item: item[1])
    n = len(scores)
    while n > 0:
    	print(scores[n-1])
    	n = n-1
    print()
    
    # MY TEST FUNCTIONS
    print("MY TEST 1")
    test1_seq, test1_header = read_fasta("./data/test_seq2.fa")
    test2_seq, test2_header = read_fasta("./data/test_seq1.fa")    
    alignment_score, seq1_alignment, seq2_alignment = test_NW.align(test1_seq, test2_seq)
    print(alignment_score)
    print(seq1_alignment)
    print(seq2_alignment)
    print()
    
    print("MY TEST 2")
    test3_seq, test1_header = read_fasta("./data/test_seq3.fa")
    test4_seq, test2_header = read_fasta("./data/test_seq4.fa")
    alignment_score, seq1_alignment, seq2_alignment = test_NW.align(test3_seq, test4_seq)
    print(alignment_score)
    print(seq1_alignment)
    print(seq2_alignment)
    

if __name__ == "__main__":
    main()
