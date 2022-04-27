#=
Set of functions to reconstruct vectors and matrices using mapping matrices, A and B
=#
using SparseArrays
using LinearAlgebra

#=
reconstruct_vec takes a reduced vector, v_reduced, and inserts the elements in v_subtracted in the appropriate locations to reconstruct the complete vector, v
    e.g. if v = { 5    and   inds = [2,4,5],   v_reduced = { 6    and   v_subtracted = { 5
                  6                                          8                           7 }
                  7                                          9 }
                  8
                  9 }

                  (inds_subtracted = [1,3])

    Therefore reconstruc_vec(v_reduced,v_subtracted,inds) will return v by performing the following linear algebraic operation:
        v = A × v_reduced + B × v_subtracted,
        where A = [ 0  0  0    and   B = [ 1  0
                    1  0  0                0  0
                    0  0  0                0  1
                    0  1  0                0  0
                    0  0  1 ]              0  0 ]
            
            (A and B are created by determining the position of the 1s in each matrix and calling a sparse function)

        Therefore:
            [ 0  0  0             [ 1  0                { 0     { 5       { 5
              1  0  0     { 6       0  0     { 5          6       0         6
              0  0  0   ×   8   +   0  1   ×   7 }   =    0   +   7    =    7
              0  1  0       9 }     0  0                  8       0         8
              0  0  1 ]             0  0 ]                9 }     0 }       9 }

Nomenclature:
    v = the complete vector, where the elements of v_reduced and v_subtracted are correctly reassembled
    v_reduced = reduced vector
    v_subtracted = the vector of elements that were removed(subtracted) from v to create v_reduced
    inds = list of element indices in v that comprise v_reduced
=#
function reconstruct_vec(v_reduced,v_subtracted,inds)
    len_v = length(v_reduced)+length(v_subtracted)                  # Length of v = length of v_reduced + length of v_subtracted
    inds_reduced = sort(inds)                                       # Sort the list of element indices in numerical order. This determines the i-position of the 1s that appear in the A matrix
    inds_subtracted = setdiff(collect(1:len_v),inds_reduced)        # Determine the indices of subtracted elements. This determines the i-position of the 1s that appear in the B matrix
    len_reduced = length(inds_reduced)                              # Length of the reduced vector
    len_subtracted = length(inds_subtracted)                        # Length of the subtracted vector
    jA = collect(1:len_reduced)                                     # The j-position of the 1s that appear in the A matrix are simply in numerical order
    sA = ones(len_reduced)                                          # Entries of the A matrix are all 1s
    A = sparse(inds_reduced,jA,sA,len_v,len_reduced)                # Create the A matrix
    jB = collect(1:len_subtracted)                                  # The j-position of the 1s that appear in the B matrix are simply in numerical order
    sB = ones(len_subtracted)                                       # Entries of the B matrix are all 1s
    B = sparse(inds_subtracted,jB,sB,len_v,len_subtracted)          # Create the B matrix
    v = A*v_reduced + B*v_subtracted                                # Perform the linear algebraic operation to reconstruct v
    return v
end

#=
reconstruct_mat takes a reduced matrix, m_reduced, and reallocates its elements to the correct position and inserts m_subtracted to reconstruct the complete matrix, m
    e.g. if m = [ 1  6 11 16 21    and   inds = [1,3,4],   m_reduced = [ 1 11 16    and   m_subtracted = [ 0  6  0  0 21    thus m_vec_subtracted = { 2 5 6 7 8 9 10 12 15 17 20 21 22 23 24 25 }^T = 16 × 1 vector
                  2  7 12 17 22                                          3 13 18                           2  7 12 17 22    (m_vec_subtracted_mod = { 0 2 0 0 5 6 7 8 9 10 0 12 0 0 15 0 17 0 0 20 21 22 23 24 25 }^T = 25 × 1 vector)
                  3  8 13 18 23                                          4 14 19 ]                         0  8  0  0 23
                  4  9 14 19 24                                                                            0  9  0  0 24
                  5 10 15 20 25 ]                                                                          5 10 15 20 25 ]

                  (inds_subtracted = [2,5])

    Therefore reconstruct_mat(m_reduced,m_subtracted,inds) will return m by performing the following linear algebraic operation:
        m = A × m_reduced × A^T + m_subtracted,
            where m_subtracted = reshape(m_vec_subtracted_mod,5,5),
                and m_vec_subtracted_mod = B × m_vec_subtracted = 25 × 1 vector
        where A = [ 1  0  0    and   B = 25 × 25 matrix (too large to list here, however it functions in the same manner as for reconstruct_vec by insterting 0s in the appropriate locations)
                    0  0  0
                    0  1  0
                    0  0  1
                    0  0  0 ]

            (A and B are created by determining the position of the 1s in each matrix and calling a sparse function)

        Therefore:
            [ 1  0  0
              0  0  0     [ 1 11 16     [ 1  0  0  0  0
              0  1  0   ×   3 13 18   ×   0  0  1  0  0   + reshape(B × m_vec_subtracted,5,5)
              0  0  1       4 14 19 ]     0  0  0  1  0 ]
              0  0  0 ]

          = [ 1  0 11 16  0     [ 0  6  0  0 21     [ 1  6 11 16 21
              0  0  0  0  0       2  7 12 17 22       2  7 12 17 22
              3  0 13 18  0   +   0  8  0  0 23   =   3  8 13 18 23
              4  0 14 19  0       0  9  0  0 24       4  9 14 19 24
              0  0  0  0  0 ]     5 10 15 20 25 ]     5 10 15 20 25 ]

Nomenclature:
    m = the complete matrix, where the elements of m_reduced and m_subtracted are correctly reassembled
    m_reduced = reduced matrix
    m_subtracted = matrix of elements removed(subtracted) from m to create m_reduced
    inds = list of row/column indices that selects m_reduced
=#
function reconstruct_mat(m_reduced,m_vec_subtracted,inds)
    len_m = Int(sqrt(length(m_reduced)+length(m_vec_subtracted)))   # Length of rows/columns of m = length of rows/columns of m_reduced + length of m_vec_subtracted
    inds_reduced = sort(inds)                                       # Sort the list of row/column indices for m_reuced in numerical order. This determines the i-position of the 1s that appear in the A matrix
    inds_subtracted = setdiff(collect(1:len_m),inds_reduced)        # Determine the row/column indices of subtracted elements. This determines the i-position of the 1s that appear in the B matrix
    len_reduced = length(inds_reduced)                              # Length of rows/columns of the reduced matrix
    jA = collect(1:len_reduced)                                     # The j-position of the 1s that appear in the A matrix are simply in numerical order
    sA = ones(len_reduced)                                          # Entries of the A matrix are all 1s
    A = sparse(inds_reduced,jA,sA,len_m,len_reduced)                # Create the A matrix
    iB = []                                                         # Initialise the list of i-positions of 1s that appear in B
    for ind in inds_subtracted                                      # For every row/column index, inds, of subtracted elements
        iB = union(iB,ind:len_m:len_m^2,(ind-1)*len_m+1:ind*len_m)  #= Add inds+n*len_m (for all natural numbers n) and (inds-1)*len_m+k*1 (for k = 1,2,3,..,len_m)
    end                                                                 This creates a list of all element numbers of the non-zero elements in m_subtracted =#
    iB = sort(iB)                                                   # Sort the list
    jB = collect(1:length(iB))                                      # The j-position of the 1s that appear in the A matrix are simply in numerical order
    sB = ones(length(iB))                                           # Entries of the A matrix are all 1s
    B = sparse(iB,jB,sB)                                            # Create the A matrix
    m = A*m_reduced*A'+reshape(B*m_vec_subtracted,len_m,len_m)      # Perform the linear algebraic operation to reconstruct m
    return m
end

# Author: Edward Street
# Date: 27/04/2022
