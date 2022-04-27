#=
Mask functions that return a binary mask matrix and the list of controlled elements for different shapes
Nomenclature:
    nelx = no. elements in the x-direction
    nely = no. elements in the y-direction
    h = x-position of the centre point as a ratio along the x-dimension
    k = y-position of the centre point as a ratio along the y-dimension
    a = length of horizontal semi-axis as a ratio of nelx
    b = length of vertical semi-axis as a ratio of nely
    r = length of radius as a ratio of nely
    l = length from centre to vertical edge as a ratio of nelx
    w = length from centre to horizontal edge as a ratio of nely
    ar = aspect ratio of nelx with respect to nely
=#
function ellipse(nelx,nely,h,k,a,b)
    mask = zeros(Int,nely,nelx)                 # Initialise mask matrix with zeros
    ctr = 0                                     # Initialise element number counter
    list = []                                   # Initialise list of controlled element numbers
    sdf(x,y) = ((x-h)/a)^2+((y-k)/b)^2-1        # sdf for an ellipse
    for j in 1:nelx                             # For every column
        for i in 1:nely                         # For every row (thus for every element, following element numbering convention)
            ctr += 1                            # Increment element number
            x = (j-0.5)/nelx                    # Shape function to map j ∈ [1,nelx] → x ∈ [0,1]
            y = (i-0.5)/nely                    # Shape function to map i ∈ [1,nely] → y ∈ [0,1]
            if sdf(x,y) <= 0                    # If sdf is negative (inside the ellipse) or 0 (on the border)
                mask[i,j] = 1                   # Update element in mask matrix with 1
                list = cat(list, ctr, dims=1)   # Append element number to the list
            end
        end
    end
    return mask, list
end

function circle(nelx,nely,h,k,r)
    mask = zeros(Int,nely,nelx)                 # Initialise mask matrix with zeros
    ctr = 0                                     # Initialise element number counter
    list = []                                   # Initialise list of controlled element numbers
    ar = nelx/nely                              # Aspect ratio to convert lengths as a ratio of nelx to lengths as a ratio of nely
    sdf(x,y) = (x-ar*h)^2+(y-k)^2-r^2           # sdf for a circle, h is converted back to a length as a ratio of nelx by multiplying the aspect ratio
    for j in 1:nelx                             # For every column
        for i in 1:nely                         # For every row (thus for every element, following element numbering convention)
            ctr += 1                            # Increment element number
            x = (j-0.5)/nely                    # Shape function to map j ∈ [1,nely] → x ∈ [0,1], there fore h and r are calculated intermediately as values as a ratio of nely
            y = (i-0.5)/nely                    # Shape function to map i ∈ [1,nely] → y ∈ [0,1]
            if sdf(x,y) <= 0                    # If sdf is negative (inside the ellipse) or 0 (on the border)
                mask[i,j] = 1                   # Update element in mask matrix with 1
                list = cat(list, ctr, dims=1)   # Append element number to the list
            end
        end
    end
    return mask, list
end

function rectangle(nelx,nely,h,k,l,w)
    mask = zeros(Int,nely,nelx)                 # Initialise mask matrix with zeros
    ctr = 0                                     # Initialise element number counter
    list = []                                   # Initialise list of controlled element numbers
    sdf(x,y) = sqrt((max(0,abs(x-h)-l))^2+(max(0,abs(y-k)-w))^2) # sdf for a rectangle (Inigo Quilez's formulation)
    for j in 1:nelx                             # For every column
        for i in 1:nely                         # For every row (thus for every element, following element numbering convention)
            ctr += 1                            # Increment element number
            x = (j-0.5)/nelx                    # Shape function to map j ∈ [1,nelx] → x ∈ [0,1]
            y = (i-0.5)/nely                    # Shape function to map i ∈ [1,nely] → y ∈ [0,1]
            if sdf(x,y) <= 0                    # If sdf is negative (inside the ellipse) or 0 (on the border)
                mask[i,j] = 1                   # Update element in mask matrix with 1
                list = cat(list, ctr, dims=1)   # Append element number to the list
            end
        end
    end
    return mask, list
end

function square(nelx,nely,h,k,w)
    mask = zeros(Int,nely,nelx)                 # Initialise mask matrix with zeros
    ctr = 0                                     # Initialise counter for element number
    list = []                                   # Initialise list of controlled element numbers
    ar = nelx/nely                              # Aspect ratio to convert lengths as a ratio of nelx to lengths as a ratio of nely
    sdf(x,y) = sqrt((max(0,abs(x-ar*h)-w))^2+(max(0,abs(y-k)-w))^2) # sdf for a square (Inigo Quilez's formulation), h is converted back to a length as a ratio of nelx by multiplying the aspect ratio
    for j in 1:nelx                             # For every column
        for i in 1:nely                         # For every row (thus for every element, following element numbering convention)
            ctr += 1                            # Increment element number
            x = (j-0.5)/nely                    # Shape function to map j ∈ [1,nely] → x ∈ [0,1], there fore h and r are calculated intermediately as values as a ratio of nely
            y = (i-0.5)/nely                    # Shape function to map i ∈ [1,nely] → y ∈ [0,1]
            if sdf(x,y) <= 0                    # If sdf is negative (inside the ellipse) or 0 (on the border)
                mask[i,j] = 1                   # Update element in mask matrix with 1
                list = cat(list, ctr, dims=1)   # Append element number to the list
            end
        end
    end
    return mask, list
end

#= 
meshoverwrite is used to overwrite elements of a matrix, a, with the nonzero elements of another matrix, b, with the same dimensions
    e.g.
    if a = [ 1  1  1  0    and   b = [ 0  0  0  0    therefore mb = [ 0  0  0  0    then meshoverwrite will return the matrix: [ 1  1  1  0
             1  1  1  0                0  0  2  2                     0  0  1  1                                                 1  1  2  2
             1  1  1  0                0  0  2  2                     0  0  1  1                                                 1  1  2  2
             0  0  0  0 ]              0  0  2  2 ]                   0  0  1  1 ]                                               0  0  2  2 ]

Nomenclature:
    a = matrix a
    b = matrix b
    mb = binary mask matrix for b, which has a 1 if the corresponding element in b is non-zero, otherwise has a 0
=#
function meshoverwrite(a,b,mb) # Overwrites elements of a matrix, a, with the nonzero elements of another matrix, b
    fa = float(a)                               # Convert all elements of matrix a into floats
    fb = float(b)                               # Convert all elements of matrix b into floats
    fa = fa.*iszero.(mb)+fb                     # First change the elements in a which will be overwritten into 0, and then add b to fill in those spots
    return fa
end

# Author: Edward Street
# Date: 27/04/2022
