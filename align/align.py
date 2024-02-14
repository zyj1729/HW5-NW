# Importing Dependencies
import numpy as np
from typing import Tuple

# Defining class for Needleman-Wunsch Algorithm for Global pairwise alignment
class NeedlemanWunsch:
    """ Class for NeedlemanWunsch Alignment

    Parameters:
        sub_matrix_file: str
            Path/filename of substitution matrix
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty

    Attributes:
        seqA_align: str
            seqA alignment
        seqB_align: str
            seqB alignment
        alignment_score: float
            Score of alignment from algorithm
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty
    """
    def __init__(self, sub_matrix_file: str, gap_open: float, gap_extend: float):
        # Init alignment and gap matrices
        self._align_matrix = None
        self._gapA_matrix = None
        self._gapB_matrix = None

        # Init matrices for backtrace procedure
        self._back = None
        self._back_A = None
        self._back_B = None

        # Init alignment_score
        self.alignment_score = 0

        # Init empty alignment attributes
        self.seqA_align = ""
        self.seqB_align = ""

        # Init empty sequences
        self._seqA = ""
        self._seqB = ""

        # Setting gap open and gap extension penalties
        self.gap_open = gap_open
        assert gap_open < 0, "Gap opening penalty must be negative."
        self.gap_extend = gap_extend
        assert gap_extend < 0, "Gap extension penalty must be negative."

        # Generating substitution matrix
        self.sub_dict = self._read_sub_matrix(sub_matrix_file) # substitution dictionary

    def _read_sub_matrix(self, sub_matrix_file):
        """
        DO NOT MODIFY THIS METHOD! IT IS ALREADY COMPLETE!

        This function reads in a scoring matrix from any matrix like file.
        Where there is a line of the residues followed by substitution matrix.
        This file also saves the alphabet list attribute.

        Parameters:
            sub_matrix_file: str
                Name (and associated path if not in current working directory)
                of the matrix file that contains the scoring matrix.

        Returns:
            dict_sub: dict
                Substitution matrix dictionary with tuple of the two residues as
                the key and score as value e.g. {('A', 'A'): 4} or {('A', 'D'): -8}
        """
        with open(sub_matrix_file, 'r') as f:
            dict_sub = {}  # Dictionary for storing scores from sub matrix
            residue_list = []  # For storing residue list
            start = False  # trigger for reading in score values
            res_2 = 0  # used for generating substitution matrix
            # reading file line by line
            for line_num, line in enumerate(f):
                # Reading in residue list
                if '#' not in line.strip() and start is False:
                    residue_list = [k for k in line.strip().upper().split(' ') if k != '']
                    start = True
                # Generating substitution scoring dictionary
                elif start is True and res_2 < len(residue_list):
                    line = [k for k in line.strip().split(' ') if k != '']
                    # reading in line by line to create substitution dictionary
                    assert len(residue_list) == len(line), "Score line should be same length as residue list"
                    for res_1 in range(len(line)):
                        dict_sub[(residue_list[res_1], residue_list[res_2])] = float(line[res_1])
                    res_2 += 1
                elif start is True and res_2 == len(residue_list):
                    break
        return dict_sub

    def align(self, seqA: str, seqB: str) -> Tuple[float, str, str]:
        """
        This function performs global sequence alignment of two sequences using the Needleman-Wunsch algorithm.
        It initializes matrices to keep track of alignment scores, backtracking directions, and gap penalties, then
        fills these matrices based on the scoring scheme. Finally, it calls the _backtrace method to construct the
        alignment from the matrices.

        Parameters:
            seqA (str): The first sequence to be aligned.
            seqB (str): The second sequence to be aligned.

        Returns:
            Tuple[float, str, str]: A tuple containing the alignment score and the aligned sequences.
        """

        # Preparation: reset alignment and score for potential re-use of this method
        self.seqA_align = ""
        self.seqB_align = ""
        self.alignment_score = 0

        # Keep original sequences for backtracing
        self._seqA = seqA
        self._seqB = seqB

        # Initialize matrices for scores and backtracking
        # Score matrix with gap penalties for first row/column
        tmp_ali_scores = np.zeros([len(seqA) + 1, len(seqB) + 1])
        tmp_ali_scores[0] = [self.gap_open + i * self.gap_extend if i else 0 for i in range(len(seqB) + 1)]
        for i in range(1, len(seqA) + 1):
            tmp_ali_scores[i][0] = self.gap_open + i * self.gap_extend
        self._ali_scores = tmp_ali_scores

        # Backtracking matrix to store directions ("<" for left, "^" for up, "`" for diagonal)
        tmp_bt = np.zeros([len(seqA) + 1, len(seqB) + 1], dtype="str")
        tmp_bt[0] = ["*"] + ["<" for i in range(1, len(seqB) + 1)]
        for i in range(1, len(seqA) + 1):
            tmp_bt[i][0] = "^"
        self._bt = tmp_bt

        # Gap matrix to indicate whether a gap was opened
        tmp_gap = np.zeros([len(seqA) + 1, len(seqB) + 1], dtype="str")
        tmp_gap[0] = ["*"] + ["g" for i in range(1, len(seqB) + 1)]
        for i in range(1, len(seqA) + 1):
            tmp_gap[i][0] = "g"
        self._gap = tmp_gap

        # Main loop to fill in the matrices based on dynamic programming
        for i in range(1, len(seqA) + 1):
            for j in range(1, len(seqB) + 1):
                # Calculate scores for moving from top, left, or diagonally (match/mismatch)
                sTop = self._ali_scores[i - 1, j] + (self.gap_open + self.gap_extend if self._gap[i - 1, j] != "g" else self.gap_extend)
                sLeft = self._ali_scores[i, j - 1] + (self.gap_open + self.gap_extend if self._gap[i, j - 1] != "g" else self.gap_extend)
                sMatch = self._ali_scores[i - 1, j - 1] + self.sub_dict[(seqA[i - 1], seqB[j - 1])]

                # Choose the highest score and update matrices accordingly
                amax = np.argmax([sTop, sLeft, sMatch])
                if amax == 0:
                    self._ali_scores[i, j] = sTop
                    self._bt[i, j] = "^"
                    self._gap[i, j] = "g"
                elif amax == 1:
                    self._ali_scores[i, j] = sLeft
                    self._bt[i, j] = "<"
                    self._gap[i, j] = "g"
                elif amax == 2:
                    self._ali_scores[i, j] = sMatch
                    self._bt[i, j] = "`"
                    self._gap[i, j] = "*"

        # After filling the matrices, perform backtracing to construct the aligned sequences
        return self._backtrace()

    def _backtrace(self) -> Tuple[float, str, str]:
        """
        This private method backtraces through the backtracking matrix to construct the aligned sequences
        and calculates the final alignment score. It starts from the bottom-right corner of the matrix and
        moves according to the backtracking directions until it reaches the top-left corner.

        Parameters:
            None

        Returns:
            Tuple[float, str, str]: A tuple containing the alignment score and the aligned sequences.
        """

        # Initialize variables for aligned sequences
        A_align = ""
        B_align = ""
        m = len(self._seqA)
        n = len(self._seqB)

        # Backtrace from bottom-right to top-left
        while self._bt[m, n] != "*":
            if self._bt[m, n] == "^":  # Move up
                A_align += self._seqA[m - 1]
                B_align += "-"
                m -= 1
            elif self._bt[m, n] == "<":  # Move left
                A_align += "-"
                B_align += self._seqB[n - 1]
                n -= 1
            elif self._bt[m, n] == "`":  # Move diagonally
                A_align += self._seqA[m - 1]
                B_align += self._seqB[n - 1]
                m -= 1
                n -= 1

        # Reverse the aligned sequences as the backtracing starts from the end
        self.seqA_align = A_align[::-1]
        self.seqB_align = B_align[::-1]

        # Final alignment score is in the bottom-right corner of the score matrix
        self.alignment_score = self._ali_scores[-1, -1]

        return (self.alignment_score, self.seqA_align, self.seqB_align)

def read_fasta(fasta_file: str) -> Tuple[str, str]:
    """
    DO NOT MODIFY THIS FUNCTION! IT IS ALREADY COMPLETE!

    This function reads in a FASTA file and returns the associated
    string of characters (residues or nucleotides) and the header.
    This function assumes a single protein or nucleotide sequence
    per fasta file and will only read in the first sequence in the
    file if multiple are provided.

    Parameters:
        fasta_file: str
            name (and associated path if not in current working directory)
            of the Fasta file.

    Returns:
        seq: str
            String of characters from FASTA file
        header: str
            Fasta header
    """
    assert fasta_file.endswith(".fa"), "Fasta file must be a fasta file with the suffix .fa"
    with open(fasta_file) as f:
        seq = ""  # initializing sequence
        first_header = True
        for line in f:
            is_header = line.strip().startswith(">")
            # Reading in the first header
            if is_header and first_header:
                header = line.strip()  # reading in fasta header
                first_header = False
            # Reading in the sequence line by line
            elif not is_header:
                seq += line.strip().upper()  # generating full sequence
            # Breaking if more than one header is provided in the fasta file
            elif is_header and not first_header:
                break
    return seq, header
