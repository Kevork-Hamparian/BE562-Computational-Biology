import sys
import numpy as np

#   dictionary for base to index
base_idx = {'A': 0, 'G': 1, 'C': 2, 'T': 3 }

#   
PTR_NONE, PTR_GAP1, PTR_GAP2, PTR_BASE = 0, 1, 2, 3


def seqalignDP(seq1, seq2, subst_matrix, gap_penalty):
    """
    Return the score of the optimal Needdleman-Wunsch alignment for seq1
    and seq2.
    Note: gap_penalty should be positive (it is subtracted)
    """
    
#    optional simple test sequences to make sure algorithm works
#    seq1='AGCTG'
#    seq2='TACGCAG'
    
    
    
    F = [[0 for j in range(len(seq2)+1)] for i in range(len(seq1)+1)]
    TB = [[PTR_NONE for j in range(len(seq2)+1)] for i in range(len(seq1)+1)]

    # initialize dynamic programming table for Needleman-Wunsch alignment
    # (Durbin p.20)
    
   
    
    # prints initial gene sequences
    print(seq1)
    print(end='\n')
    print(seq2)
    print(end='\n')
    
###         Part a  
    ''' to be un-commented if you are trying to find the final 
        score of the sequence alignment'''
        
#    for i in range(1, len(seq1)+1):
#        F[i][0] = 0 - i*gap_penalty
#        TB[i][0] = PTR_GAP2  # indicates a gap in seq2
#        
#    for j in range(1, len(seq2)+1):
#        F[0][j] = 0 - j*gap_penalty
#        TB[0][j] = PTR_GAP1  # indicates a gap in seq1
#        
#    for i in range(1, len(seq1)+1):
#        for j in range(1, len(seq2)+1):
#            
#
#            gap2=F[i-1][j]-gap_penalty
#            gap1=F[i][j-1]-gap_penalty
#            
#            l1=base_idx[seq1[i-1]]
#            l2=base_idx[seq2[j-1]]
#            match_coef=subst_matrix[l2][l1]
#            
#            diag=F[i-1][j-1]+match_coef
#            
#            F[i][j]=max(gap1,gap2,diag)
#            if      F[i][j]==gap1:
#                       TB[i][j]=PTR_GAP1
#            elif    F[i][j]==gap2:
#                       TB[i][j]=PTR_GAP2
#            elif    F[i][j]==diag:
#                       TB[i][j]=PTR_BASE
#            else:
#                       TB[i][j]=PTR_NONE
                       
##             Part B   
    """
        to be un-commented if you are trying to find the distance
        between the two sequences
    """       
    for i in range(1, len(seq1)+1):
        F[i][0] = 0 + i*gap_penalty
        TB[i][0] = PTR_GAP2  # indicates a gap in seq2
        
    for j in range(1, len(seq2)+1):
        F[0][j] = 0 + j*gap_penalty
        TB[0][j] = PTR_GAP1  # indicates a gap in seq1
        
    for i in range(1, len(seq1)+1):
        for j in range(1, len(seq2)+1):
            gap2=F[i-1][j]+gap_penalty
            gap1=F[i][j-1]+gap_penalty
            
            l1=base_idx[seq1[i-1]]
            l2=base_idx[seq2[j-1]]
            match_coef=subst_matrix[l2][l1]
            
            diag=F[i-1][j-1]+match_coef
            
            F[i][j]=min(gap1,gap2,diag)
            
            if      F[i][j]==gap1:
                        TB[i][j]=PTR_GAP1
            elif    F[i][j]==gap2:
                        TB[i][j]=PTR_GAP2
            elif    F[i][j]==diag:
                        TB[i][j]=PTR_BASE

    #print(F)
    

    # YOUR CODE HERE
    # Fill in the dynamic programming tables F and TB, starting at [1][1]
    # Hints: The first row and first column of the table F[i][0] and F[0][j]
    # contain dummy values.
    #  (see for illustration Durbin p.21, Figure 2.5, but be careful what you
    #   think of as rows and what you think of as columns)
    #  Hence, the bases corresponding to F[i][j] are actually seq1[i-1] and
    #  seq2[j-1].
    

    #  Use the dictionary base_idx to convert from the character to an index to
    #   look up entries of the substitution matrix.
    #  To get started, you can complete and run the algorithm filling in only
    #    F, and then figure out how to do TB.
    
#    for i in range(len(F)):
#        print(F[i],'\n')   
#        
#    print(TB)

    return F[len(seq1)][len(seq2)], F, TB


def traceback(seq1, seq2, TB):
#   Create's traceback matrix where you can follow the optimal path    
    '''
        seq1 and seq2 are to be uncommented only to 
        test simple sequences
    '''
#    seq1='AGCTG'
#    seq2='TACGCAG'
    
    s1 = ""
    s2 = ""

    i = len(seq1)
    j = len(seq2)
    

    
    while TB[i][j] != PTR_NONE:
        if TB[i][j] == PTR_BASE:
            s1 = seq1[i-1] + s1
            s2 = seq2[j-1] + s2
            i = i - 1
            j = j - 1
        elif TB[i][j] == PTR_GAP1:
            s1 = '-' + s1
            s2 = seq2[j-1] + s2
            j = j - 1
        elif TB[i][j] == PTR_GAP2:
            s1 = seq1[i-1] + s1
            s2 = '-' + s2
            i = i - 1
        else:
            assert False
    
    print(TB)
    return s1, s2


def readSeq(filename):
    """Reads in a FASTA sequence. Assumes one sequence in the file"""
    seq = []

    with open(filename, "r") as f:
        for line in f:
            if line.startswith(">"):
                continue
            seq.append(line.rstrip().upper())

    return "".join(seq)

# Substituation matrix and gap_penalty
'''
    to be uncommented if you would like to score the optimal sequence
'''
#S = [
#    # A   G   C   T
#    [ 3, -1, -2, -2],  # A
#    [-1,  3, -2, -2],  # G
#    [-2, -2,  3, -1],  # C
#    [-2, -2, -1,  3]   # T
#]

# Substitution matrix for distance calc
S = [
    # A   G   C   T
    [ 0,  1,  2,  2],  # A
    [ 1,  0,  2,  2],  # G
    [ 2,  2,  0,  1],  # C
    [ 2,  2,  1,  0]   # T
]

#   The penalty of inserting a gap into one of the sequences
gap_penalty = 4


def main():
    # parse command line
    if len(sys.argv) < 3:
        print("Usage: {0} <FASTA 1> <FASTA 2>".format(sys.argv[0]))
        sys.exit(1)

    file1 = sys.argv[1]
    file2 = sys.argv[2]

    seq1 = readSeq(file1)
    seq2 = readSeq(file2)

    score, F, TB = seqalignDP(seq1, seq2, S, gap_penalty)

    print("Score: {0}".format(score))

    s1, s2 = traceback(seq1, seq2, TB)
    print(s1)
    print(s2)

if __name__ == "__main__":
    main()
