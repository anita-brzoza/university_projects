
!!! Unfinished

This folder contain the [Smith-Waterman_algorithm.py](/university_projects/Smith-Waterman_algorithm/Smith-Waterman_algorithm.py) file with a function using the Smith-Waterman algorithm for local matching two sequences.

Smith-Waterman algorithm - Needelman-Wunsh algorithm in local matching version:

                   H(i - 1, j - 1) + S(a_i, b_j)   
    H(i,j) = max { H(i - 1, j) - g
                   H(i, j - 1) - g
                   0

    where:
        S(a, b) - similarity matrix beetween sequence A and B
        H - scoring matrix
        g - responsible for the penalty function