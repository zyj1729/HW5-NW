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
    NW = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)
    score_gg, hs_align, gg_align = NW.align(hs_seq, gg_seq)
    score_mm, _, mm_align = NW.align(hs_seq, mm_seq)
    score_br, _, br_align = NW.align(hs_seq, br_seq)
    score_tt, _, tt_align = NW.align(hs_seq, tt_seq)
    species = ["Gallus_gallus", "Mus_musculus", "Balaeniceps_rex", "tursiops_truncatus"]
    scores = [score_gg, score_mm, score_br, score_tt]
    rank = sorted([(scores[i], species[i]) for i in range(len(species))], reverse = True)
    for i in rank:
        print(i[1])

    # TODO print all of the alignment score between each species BRD2 and human BRD2
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    for i in rank:
        print("Score for human BRD2 and " + i[1] + " BRD2:", i[0])
    

if __name__ == "__main__":
    main()
