using SparseArrays
using LinearAlgebra
using Statistics
using Printf
using Plots

function top88(nelx,nely,volfrac,penal,rmin,ft)
    # Material Properties
    E0 = 1
    Emin = 1e-9
    nu = 0.3
    # Prepare Finite Element Analysis
    A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12]
    A12 = [-6 -3  0  3; -3 -6 -3  -6;  0 -3 -6  3;  3 -6  3 -6]
    B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4]
    B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2]
    KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11])
    nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx)
    edofVec = reshape(2*nodenrs[1:end-1,1:end-1].+1,nelx*nely,1)
    edofMat = repeat(edofVec,1,8)+repeat([0 1 2*nely.+[2 3 0 1] -2 -1],nelx*nely,1)
    iK = reshape(kron(edofMat,ones(Int,8,1))',64*nelx*nely,1)
    jK = reshape(kron(edofMat,ones(Int,1,8))',64*nelx*nely,1)
    # Define Loads and Supports (Half MBB-Beam)
    F = sparse([2],[1],[-1],2*(nely+1)*(nelx+1),1)
    U = zeros(2*(nely+1)*(nelx+1),1)
    fixeddofs = union(collect(1:2:2*(nely+1)),[2*(nelx+1)*(nely+1)])
    alldofs = collect(1:2*(nelx+1)*(nely+1))
    freedofs = setdiff(alldofs,fixeddofs)
    # Prepare Filter
    iH = ones(Int,nelx*nely*(2*(ceil(Int,rmin)-1)+1)^2,1)
    jH = ones(Int,size(iH))
    sH = zeros(size(iH))
    k = 0
    for i1 in 1:nelx
        for j1 in 1:nely
            e1 = (i1-1)*nely+j1
            for i2 in max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
                for j2 in max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                  e2 = (i2-1)*nely+j2
                  k = k+1
                  iH[k] = e1
                  jH[k] = e2
                  sH[k] = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2))
                end
            end
        end
    end
    H = sparse(vec(iH),vec(jH),vec(sH))
    Hs = sparse(sum(H,dims=2))
    # Initialise Iteration
    x = repeat([volfrac],nely,nelx)
    xPhys = copy(x)
    loop = 0
    change = 1
    # Start Iteration
    while change > 0.01
        loop += 1
        # FE-Analysis
        sK = reshape(KE[:]*(Emin.+xPhys[:]'.^penal*(E0-Emin)),64*nelx*nely,1)
        K = sparse(vec(iK),vec(jK),vec(sK)); K = (K+K')/2
        U[freedofs] = cholesky(K[freedofs, freedofs])\Vector(F[freedofs])
        # Objective Function and Sensitivity Analysis
        ce = reshape(sum((U[edofMat]*KE).*U[edofMat],dims=2),nely,nelx)
        c = sum(sum((Emin.+xPhys.^penal*(E0-Emin)).*ce))
        dc = -penal*(E0-Emin)*xPhys.^(penal-1).*ce
        dv = ones(nely,nelx)
        # Filtering/Modification of Sensitivities
        if ft == 1
            dc[:] = H*(x[:].*dc[:])./Hs./max.(fill(1e-3,size(x[:])),x[:])
        elseif ft == 2
            dc[:] = H*(dc[:]./Hs)
            dv[:] = H*(dv[:]./Hs)
        end
        # Optimality Criteria Update of Design Variables and Physical Densities
        l1 = 0; l2 = 1e9; move = 0.2; xnew = zeros(nely,nelx)
        while (l2-l1)/(l1+l2) > 1e-3
            lmid = 0.5*(l2+l1)
            xnew = max.(fill(0,nely,nelx),max.(x.-move,min.(fill(1,nely,nelx),min.(x.+move,x.*sqrt.(-dc./dv/lmid)))))
            if ft == 1
                xPhys = xnew
            elseif ft == 2
                xPhys[:] = (H*xnew[:])./Hs
            end
            if sum(xPhys[:]) > volfrac*nelx*nely; l1 = lmid; else l2 = lmid; end
        end
        change = maximum(abs.(xnew[:]-x[:]))
        x = xnew
        # Print Results
        @printf(" It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n",loop,c,
            mean(xPhys[:]),change)
        # Plot Densities
        display(heatmap(xPhys, xlim=(0,nelx+1), ylim=(0,nely+1), color=:Greys, clims=(0,1), aspect_ratio=1, yflip=true))
    end
end

# Author: Edward Street
# Date: 27/04/2022