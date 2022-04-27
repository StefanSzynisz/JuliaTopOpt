using SparseArrays
using LinearAlgebra
using Statistics
using Printf
using Plots

function ellipse(nelx,nely,h,k,a,b)
    mask = zeros(Int,nely,nelx)
    ctr = 0
    list = []
    sdf(x,y) = ((x-h)/a)^2+((y-k)/b)^2-1
    for j in 1:nelx
        for i in 1:nely
            ctr += 1
            x = (j-0.5)/nelx
            y = (i-0.5)/nely
            if sdf(x,y) <= 0
                mask[i,j] = 1
                list = cat(list, ctr, dims=1)
            end
        end
    end
    return mask, list
end

function circle(nelx,nely,h,k,r)
    mask = zeros(Int,nely,nelx)
    ctr = 0
    list = []
    ar = nelx/nely
    for j in 1:nelx
        for i in 1:nely
            ctr += 1
            x = (j-0.5)/nely
            y = (i-0.5)/nely
            if sdf(x,y) <= 0
                mask[i,j] = 1
                list = cat(list, ctr, dims=1)
            end
        end
    end
    return mask, list
end

function rectangle(nelx,nely,h,k,l,w)
    mask = zeros(Int,nely,nelx)
    ctr = 0
    list = []
    sdf(x,y) = sqrt((max(0,abs(x-h)-l))^2+(max(0,abs(y-k)-w))^2)
    for j in 1:nelx
        for i in 1:nely
            ctr += 1
            x = (j-0.5)/nelx
            y = (i-0.5)/nely
            if sdf(x,y) <= 0
                mask[i,j] = 1
                list = cat(list, ctr, dims=1)
            end
        end
    end
    return mask, list
end

function square(nelx,nely,h,k,w)
    mask = zeros(Int,nely,nelx)
    ctr = 0
    list = []
    ar = nelx/nely
    sdf(x,y) = sqrt((max(0,abs(x-ar*h)-w))^2+(max(0,abs(y-k)-w))^2)
    for j in 1:nelx
        for i in 1:nely
            ctr += 1
            x = (j-0.5)/nely
            y = (i-0.5)/nely
            if sdf(x,y) <= 0
                mask[i,j] = 1
                list = cat(list, ctr, dims=1)
            end
        end
    end
    return mask, list
end

function meshoverwrite(a,b,mb)
    fa = float(a)
    fb = float(b)
    fa = fa.*iszero.(mb)+fb
    return fa
end

function reconstruct_vec(v_red,v_sub,inds)
    len_v = length(v_red)+length(v_sub)
    inds_reduced = sort(inds)
    inds_subtracted = setdiff(collect(1:len_v),inds_reduced)
    len_reduced = length(inds_reduced)
    len_subtracted = length(inds_subtracted)
    jA_reduced = collect(1:len_reduced)
    sA_reduced = ones(len_reduced)
    A_reduced = sparse(inds_reduced,jA_reduced,sA_reduced,len_v,len_reduced)
    jA_subtracted = collect(1:len_subtracted)
    sA_subtracted = ones(len_subtracted)
    A_subtracted = sparse(inds_subtracted,jA_subtracted,sA_subtracted,len_v,len_subtracted)
    v = A_reduced*v_red + A_subtracted*v_sub
    return v
end

function reconstruct_mat(m_reduced,m_vec_subtracted,inds)
    len_m = Int(sqrt(length(m_reduced)+length(m_vec_subtracted)))
    inds_reduced = sort(inds)
    inds_subtracted = setdiff(collect(1:len_m),inds_reduced)
    len_reduced = length(inds_reduced)
    jA_reduced = collect(1:len_reduced)
    sA_reduced = ones(len_reduced)
    A_reduced = sparse(inds_reduced,jA_reduced,sA_reduced,len_m,len_reduced)
    iB = []
    for ind in inds_subtracted
        iB = union(iB,ind:len_m:len_m^2,(ind-1)*len_m+1:ind*len_m)
    end
    iB = sort(iB)
    jB = collect(1:length(iB))
    sB = ones(length(iB))
    B = sparse(iB,jB,sB)
    m = A_reduced*m_reduced*A_reduced'+reshape(B*m_vec_subtracted,len_m,len_m)
    return m
end

function top88_dcon(nelx,nely,volfrac,penal,rmin,ft)
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
    # Call Mask Functions and Initialise Iteration
    r1, r1_list = rectangle(nelx,nely,0.75,0.5,0.05,0.3)        # Mask matrix and the list of controlled elements for a rectangle
    e1, e1_list = ellipse(nelx,nely,0.25,0.5,0.1,0.4)           # Mask matrix and the list of controlled elements for an ellipse
    r1_val = 1                                                  # Set rectangle as solid
    e1_val = 0                                                  # Set ellipse as void
    mask = r1 .| e1                                             # Merge binary mask matrices (1 for controlled/inactive elements, 0 for active elements)
    notmask = iszero.(mask)                                     # Return reverse of mask (0 for controlled/inactive elements, 1 for active elements)
    controlmesh = meshoverwrite(r1_val*r1,e1_val*e1,e1)         # Mask matrix with prescribed density values in the appropriate positions (the ellipse is overwritten on top of the rectangle)
    elements_inactive = sort(union(r1_list,e1_list))            # Array of element numbers for inactive elements
    elements_active = setdiff(reshape(1:nelx*nely,nelx*nely,1),elements_inactive) # Array of element numbers for active elements
    len_active = length(elements_active)                        # No. of active elements
    controlvol = sum(controlmesh)                               # Amount of density that is predefined by controlmesh
    new_volfrac = (volfrac*nelx*nely-controlvol)/(sum(notmask)) # Amount of remaining density to equally distribute between active elements
    x = controlmesh + new_volfrac*notmask                       # Initialise the mesh by equally redistributing new_volfrac to active elements
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
        dc_vec_reduced = dc[:][elements_active]                 # Vectorised form of dc, reduced to contain values only for active elements
        dv_vec_reduced = dv[:][elements_active]                 # Vectorised form of dv, reduced to contain values only for active elements
        x_vec_reduced = x[:][elements_active]                   # Vectorised form of x, reduced to contain values only for active elements
        x_vec_subtracted = x[:][elements_inactive]              # The vector of elements removed(subtracted) from x[:] when creating x_vec_reduced
        # Optimality Criteria Update of Design Variables and Physical Densities
        l1 = 0; l2 = 1e9; move = 0.2; xnew = zeros(nely,nelx)
        while (l2-l1)/(l1+l2) > 1e-3
            lmid = 0.5*(l2+l1)
            B = -dc_vec_reduced./dv_vec_reduced/lmid            # Calculate the Optimality Criteria updating factor using the reduced vectors
            xnew_vec_reduced = max.(zeros(len_active),max.(x_vec_reduced.-move,min.(ones(len_active),
                min.(x_vec_reduced.+move,x_vec_reduced.*sqrt.(B))))) # The vectorised form of xnew, only conataining values for active elements
            xnew_vec = reconstruct_vec(xnew_vec_reduced,x_vec_subtracted,elements_active)   # Vectorised form of xnew, reconstructed for all elements
            xnew = reshape(xnew_vec,nely,nelx)                  # Reshape xnew back to a matrix
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