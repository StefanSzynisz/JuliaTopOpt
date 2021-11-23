using SparseArrays
using LinearAlgebra

function topopt(nelx,nely)
    ndim = 2                                                                                    #2D elements
    neltot = nelx*nely                                                                          #no. total elements
    ngnodes = (nelx+1)*(nely+1)                                                                 #no. global nodes
    nlocnodes = 4                                                                               #no. local nodes (4-noded bilinear element)
    ngdofs = ndim*ngnodes                                                                       #no. global dofs
    nlocdofs = ndim*nlocnodes                                                                   #no. local dofs
    gnodemat = reshape(1:ngnodes,nely+1,nelx+1)                                                 #map of global nodes
    edofvec = reshape(2*gnodemat[1:end-1,1:end-1].+1,neltot,1)                                  #column vector of leading global node no. for each element
    edofmat = repeat(edofvec,1,nlocdofs)+repeat([0 1 2*(nely).+[2 3 0 1] -2 -1],neltot,1)       #global dofs in order for each element
    E0 = 1                                                                                      #Initial stiffness
    Emin = 1e-9                                                                                 #Minimal stiffness value for penalisation
    nu = 0.3                                                                                    #Poisson's ratio
    A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12]                                  #top-left and bottom-right in [A]
    A12 = [-6 -3  0  3; -3 -6 -3  6;  0 -3 -6  3;  3 -6  3 -6]                                  #top-right and bottom-left in [A]
    B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4]                                  #top-left and bottom-right in [B]
    B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2]                                  #top-right and bottom-left in [B]
    A = [A11 A12; A12 A11]                                                                      #[A]
    B = [B11 B12; B12 B11]                                                                      #[B]
    ke = 1/(24*(nu^2))*(A-nu*B)                                                                 #Element stiffness matrix
    K = zeros(ngdofs,ngdofs)                                                                    #Initialise global stiffness matrix with 0s
    for i in 1:neltot                                                                           #i = For every element
        for j in 1:ndim:nlocdofs                                                                #j = 1, 3, 5, 7
            for k in 1:ndim:nlocdofs                                                            #k = 1, 3, 5, 7
                K[edofmat[i,j]:edofmat[i,j+1],edofmat[i,k]:edofmat[i,k+1]] += ke[j:j+1,k:k+1]   #Fill in global stiffness matrix by assembling values in element stiffness matrix (not sparse function yet)
            end
        end
    end
    F = zeros(ngdofs,1)                                                                         #Initialise force vector
    F[2] = -1                                                                                   #Apply vertical load of -1 to node #1 (not sparse function yet)
    U = zeros(ngdofs,1)                                                                         #Initialise displacement vector
    alldofs = collect(1:ngdofs)                                                                 #Vector of all global dofs
    fixeddofs = union(collect(1:2:2*(nely+1)),[ngdofs])                                         #Vector of all global dofs that are fixed
    freedofs = setdiff(alldofs,fixeddofs)                                                       #Vector of all global dofs that are free
    return(K)
end
