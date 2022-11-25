"""
!!! Unfinished
"""

import numpy as np
from Bio import SeqIO

def reading_similarity_matrix(path_similarity_matrix):
    """Function to open and save similarity matrix"""
    try:
        with open(path_similarity_matrix, "r") as file:
            read = file.read()
            read = read.splitlines()
    except FileNotFoundError:
        print("Path to similarity matrix file is not correct or file does not exist.")
    else:
        n = len(read) - 1
        scoring_matrix_numbers = np.zeros((n,n), "i")
        scoring_matrix_alphabet = {}
        row = read[0].split()
        for j in range(0, len(row)):
            scoring_matrix_alphabet[row[j]] = j
        for i in range(1, len(row)):
            row=read[i].split()
            for j in range(1, len(row)):
                if row[j]!='':
                    scoring_matrix_numbers[i-1,j-1] = float(row[j])
                    scoring_matrix_numbers[j-1,i-1] = scoring_matrix_numbers[i-1,j-1] 
        return scoring_matrix_numbers, scoring_matrix_alphabet

def reading_sequences(path_sequence_A, path_sequence_B):
    """Function that open the fasta file with the sequence A and B"""
    try: 
        sequence_A = [str(file.seq) for file in SeqIO.parse(path_sequence_A, "fasta")]
    except FileNotFoundError:
        print("Path to sequence A fasta file is not correct or file does not exist.")
    else:
        try:
            sequence_B = [str(file.seq) for file in SeqIO.parse(path_sequence_B, "fasta")]
        except FileNotFoundError:
            print("Path to sequence B fasta file is not correct or file does not exist.")
        else:
            return sequence_A[0], sequence_B[0]

def smithwaterman_algorithm(sequence_A, sequence_B, scoring_matrix_numbers, scoring_matrix_alphabet, penalty=6):
    """
    Function that implement the Smith-Waterman algorithm - Needelman-Wunsh algorithm in local matching version:

                   H(i - 1, j - 1) + S(a_i, b_j)   
    H(i,j) = max { H(i - 1, j) - g
                   H(i, j - 1) - g
                   0

    where:
        S(a, b) - similarity matrix beetween sequence A and B
        H - scoring matrix
        g - responsible for the penalty function


    The function returns: 
        * best score of matching
        * sequence 1
        * local matching ("|" -match, ":" - mismatch, "-" - gap)
        * sequence 2
    """

    rows = len(sequence_A)+1
    columns = len(sequence_B)+1
    score = np.zeros((rows, columns))
    arrows = np.zeros((rows, columns))
        
    for i in range(1, rows):
        for j in range(1, columns):
            p1 = scoring_matrix_alphabet[sequence_A[i-1]]
            p2 = scoring_matrix_alphabet[sequence_B[j-1]]
            match = score[i-1, j-1] + scoring_matrix_numbers[p1,p2]
            vertical = score[i-1, j] - penalty
            horizontal = score[i, j-1] - penalty
            score[i,j] = max(0,vertical,horizontal,match)
            if score[i,j] == match:
                arrows[i,j] = 3
            elif score[i,j] == horizontal:
                arrows[i,j] = 2
            elif score[i,j] == vertical:
                arrows[i,j] = 1
            else:
                arrows[i,j] = 0
    
    #the best matching
    best_score = np.amax(score)
    ind_best_score = np.where(score == best_score)
    i_max = ind_best_score[0]
    j_max = ind_best_score[1]
    
    i = rows - 1
    j = columns - 1

    sequence_A = list(sequence_A)
    sequence_B = list(sequence_B)

    while i != 0 & j != 0 & i < i_max & j < j_max:
        if arrows[i, j] == 3:
            i = i - 1
            j = j - 1
        elif arrows[i, j] == 2:
            sequence_A.insert(i, '-')
            j = j - 1
        elif arrows[i, j] == 1:
            sequence_B.insert(j, '-')
            i = i - 1
        else:
            break
    
    seq_A=''
    seq_B=''
    for a in sequence_A:
        seq_A += a
    for a in sequence_B:
        seq_B += a

    local_matching = ""
    if len(sequence_A) != len(sequence_B):
        if len(sequence_A) > len(sequence_B):
            for a in range(len(sequence_B)):
                if sequence_B[a] == sequence_A[a] and sequence_A[a] != "-" and sequence_A[a] != "-":
                    local_matching_sign = "|"
                elif sequence_A[a] != "-":
                    local_matching_sign = ":"
                else:
                    local_matching_sign = " "
                local_matching += local_matching_sign
        elif len(sequence_A) < len(sequence_B):
            for a in range(len(sequence_A)):
                if sequence_A[a] == sequence_B[a] and sequence_B[a] != "-" and sequence_A[a] != "-":
                    local_matching_sign = "|"
                elif sequence_B[a] != "-":
                    local_matching_sign = ":"
                else:
                    local_matching_sign = " "
                local_matching += local_matching_sign
    else:
        for a in range(len(sequence_A)):
            if sequence_A[a] == sequence_B[a] and sequence_B[a] != "-" and sequence_A[a] != "-":
                local_matching_sign = "|"
            elif sequence_B[a] != "-":
                local_matching_sign = ":"
            else:
                local_matching_sign = " "
            local_matching += local_matching_sign

    print(f"Best score: {best_score}\n{seq_A}\n{local_matching}\n{seq_B}")


        
# ------------------------------------------------------------------------------
"""Example: """

paths = ["D:/github/university_projects/Smith-Waterman_algorithm/seqA.fasta", 
"D:/github/university_projects/Smith-Waterman_algorithm/seqB.fasta", 
"D:/github/university_projects/Smith-Waterman_algorithm/BLOSUM62.txt"]

scoring_matrix_numbers, scoring_matrix_alphabet = reading_similarity_matrix(paths[2])
sequence_A, sequence_B = reading_sequences(paths[0], paths[1])
smithwaterman_algorithm(sequence_A, sequence_B, scoring_matrix_numbers, scoring_matrix_alphabet)