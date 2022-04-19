# JuliaTopOpt
Topology Optimization using Julia Programming language

# Matrices
We use xxx for regular matrices in concept testing but sparse matrices are used for production code.

# Method of introducing tables
| Syntax      | Description |
| ----------- | ----------- |
| Header      | Title       |
| Paragraph   | Text        |

# Example of LaTeX
<img src="https://render.githubusercontent.com/render/math?math=e^{i \pi} = -1">

# top88.m Deconstruction
| Notes/Annotations |
| :---------------- |
| <ul><li>[a x b] – denoting dimension of array/matrix</li><li>dofs – degrees of freedom</li></ul> |

| Lines | Source Code | Description |
| :---: | :---------- | :---------- |
| 2     | function top88(nelx,nely,volfrac,penal,rmin,ft) | Main function call<br /> Inputs:<ul><li>nelx – no. of elements in x direction</li><li>nely – no. of elements in y direction</li><li>volfrac – volume fraction</li><li>penal – penalisation factor</li><li>rmin – filter radius</li><li>ft – filtering mode (1: Sensitivity filtering, 2: Density Filtering)</li></ul> |
| 4-6   | E0 = 1;<br /> Emin = 1e-9;<br /> nu = 0.3; | Setting material properties:<ul><li>E0 – Nominal Stiffness</li><li>Emin – Small non-zero stiffness to avoid singular matrices</li><li>nu – Poisson’s ratio</li></ul> |
| 8-11  | A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];<br /> A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];<br /> B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];<br /> B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2]; | Constructing components for element stiffness matrix for 2D 4-noded bilinear element for plane stress |
| 12    | KE = 1/(1-nu^2)/24\*([A11 A12;A12' A11]+nu\*[B11 B12;B12' B11]); | Assembly of element stiffness matrix<br />[8 x 8] |
| 13    | nodenrs = reshape(1:(1+nelx)\*(1+nely),1+nely,1+nelx); | Create global node numbers matrix (mapping the Cartesian mesh as a 2d array of global node numberings), convention: top down, starting from leftmost column<br /> [nely+1 x nelx+1] |
| 14    | edofVec = reshape(2\*nodenrs(1:end-1,1:end-1)+1,nelx\*nely,1); | Create 1D array of leading dofs for each element in edofMat (x-direction dof for the bottom left node for each element)<br /> [nelx\*nely x 1] |
| 15    | edofMat = repmat(edofVec,1,8)+repmat(\[0 1 2\*nely+\[2 3 0 1\] -2 -1\],nelx\*nely,1); | Create element dof matrix, where each row notates the global dof number for each element with convention: starting from the bottom left x-direction dof, then y-direction dof, going anticlockwise for each element, (steering matrix)<br /> [nelx\*nely x 8] |
| 16-17 | iK = reshape(kron(edofMat,ones(8,1))',64\*nelx\*nely,1);<br /> jK = reshape(kron(edofMat,ones(1,8))',64\*nelx\*nely,1); | Preparing index vectors for sparse function to construct global stiffness matrix<br /> [64\*nelx\*nely x 1] |
| 19    | F = sparse(2,1,-1,2*(nely+1)*(nelx+1),1); | Applying load of -1 at global dof 2 (therefore for top left node of elemement 1 in y-direction)<br /> [2\*(nelx+1)\*(nely+1) x 1] |
| 20    | U = zeros(2*(nely+1)*(nelx+1),1); | Initial Displacement Vector (filled with 0s)<br /> [2\*(nelx+1)\*(nely+1) x 1] |
| 21    | fixeddofs = union([1:2:2*(nely+1)\],[2*(nelx+1)*(nely+1)]); | 1D array of fixed dofs, (For MBB beam, fixed x-displacements for left side nodes due to symmetry condition and fixed y-displacement for bottom right corner node due to fixed support) |
| 22    | alldofs = [1:2*(nely+1)*(nelx+1)]; | 1D array of all global dofs<br /> [2\*(nelx+1)\*(nely+1) element array] |
| 23    | freedofs = setdiff(alldofs,fixeddofs); | 1D array of remaining free dofs |
| 25-27 | iH = ones(nelx\*nely*(2*(ceil(rmin)-1)+1)^2,1);<br /> jH = ones(size(iH));<br /> sH = zeros(size(iH)); | Initialise index vectors for creating the filter weight matrix (H, which is a collection of all values of H<sub>ei</sub>) using the sparse function<br /> [nelx\*nely x 1] |
| 28    | k = 0; | Initialise counter for filter weight matrix construction |
| 29    | for i1 = 1:nelx | For i1 = 1, 2, ... , nelx |
| 30    | &nbsp;&nbsp;&nbsp;&nbsp;for j1 = 1:nely | For j1 = 1, 2, ... , nely (therefore, sampling through every element in the mesh) |
| 31    | &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;e1 = (i1-1)*nely+j1; | Element 1 (element e) number |
| 32-33 | &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx);<br />&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely); | for all elements which the centre-to-centre distance to element e is smaller than the filter radius |
| 34    | &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;e2 = (i2-1)*nely+j2; | Element 2 (element i) number |
| 35    | &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;k = k+1; | Increment counter |
| 36    | &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;iH(k) = e1; | Place element 1 number in iH |
| 37    | &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;jH(k) = e2; | Place element 2 number in jH |
| 38    | &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2)); | Place corresponding weight value, H<sub>ei</sub> in sH |
| 43    | H = sparse(iH,jH,sH); | Create filtering weight matrix using sparse function<br /> [nelx\*nely\*(2(&lceil;rmin&rceil;-1)+1)<sup>2</sup> x nelx\*nely\*(2(&lceil;rmin&rceil;-1)+1)<sup>2</sup>] |
| 44    | Hs = sum(H,2); | Hs = column vector containing the sum of each row of the H matrix<br /> [nelx\*nely x nelx\*nely] |
| 46    | x = repmat(volfrac,nely,nelx); | Create matrix containing density values associated to each mesh cell<br /> [nely x nelx] |
| 47    | xPhys = x; | xPhys is a copy of x used to represent the physical mesh densities after filtering<br /> [nely x nelx] |
| 48    | loop = 0; | Initialise loop count |
| 49    | change = 1; | Initialise change value |
| 51    | while change > 0.01 | Begin main while loop<br /> The change value is defined as the L<sup>&infin;</sup> norm of the difference between the density matrix of the new iteration in vectorised form, xnew(:), and the density matrix of the current iteration in vectorised form, x(:), therefore the maximum absolute value within the vector components in xnew(:)-x(:) |
| 52    | loop = loop + 1; | Increment loop counter |
| 54    | sK = reshape(KE(:)\*(Emin+xPhys(:)'.^penal\*(E0-Emin)),64\*nelx\*nely,1); | Global stiffness matrix entries after modifying the element stifnesses according to the modified SIMP formulation<br /> [64\*nelx\*nely x 1] |
| 55    | K = sparse(iK,jK,sK); K = (K+K')/2; | Create the global stiffness matrix using the sparse function |
| 56    | U(freedofs) = K(freedofs,freedofs)\F(freedofs); | Solve the FEA stiffness equation for global dof displacements |
| 58    | ce = reshape(sum((U(edofMat)\*KE).*U(edofMat),2),nely,nelx); | Calculate element compliances |
| 59    | c = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce)); | Calculate the total compliance of the design |
| 60    | dc = -penal*(E0-Emin)\*xPhys.^(penal-1).*ce; | Calculate compliance sensitivities for every element<br /> [nely x nelx] |
| 61    | dv = ones(nely,nelx); | Volume sensitivities for every element (=1 by definition)<br /> [nely x nelx] |
| 63    | if ft == 1 | If the sensitivity filter is selected |
| 64    | &nbsp;&nbsp;&nbsp;&nbsp;dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:)); | Modify the compliance sensitivity according to the sensitivity filter equation |
| 65    | elseif ft == 2 | If the density filter is selected |
| 66-67    | &nbsp;&nbsp;&nbsp;&nbsp;dc(:) = H*(dc(:)./Hs);<br />&nbsp;&nbsp;&nbsp;&nbsp;dv(:) = H*(dv(:)./Hs); | Modify both the compliance and volume sensitivities according to the density filter (chain rule) equation |
| 70    | l1 = 0; l2 = 1e9; move = 0.2; | Prepare lower and upper bound values for lambda used in the bisection method, and the move limit |
| 71    | while (l2-l1)/(l1+l2) > 1e-3 | While the bisection method is not considered converged |
| 72    | &nbsp;&nbsp;&nbsp;&nbsp;lmid = 0.5*(l2+l1); | Calculate lmid, the average of l1 and l2 |
| 73    | &nbsp;&nbsp;&nbsp;&nbsp;xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(-dc./dv/lmid))))); | Compute the density matrix for the next iteration from the Optimality Criteria method updating formulation (the three case conditional statements are combined into one line of code)<br /> [nely x nelx] |
| 74    | &nbsp;&nbsp;&nbsp;&nbsp;if ft == 1 | If the sensitivity filter is selected |
| 75    | &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;xPhys = xnew; | xnew can be accepted as the physical solution, xPhys<br /> [nely x nelx] |
| 76    | &nbsp;&nbsp;&nbsp;&nbsp;elseif ft == 2 | If the density filter is selected |
| 77    | &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;xPhys(:) = (H*xnew(:))./Hs; | The densities are modified according to the density filter equation |
| 79    | &nbsp;&nbsp;&nbsp;&nbsp;if sum(xPhys(:)) > volfrac\*nelx\*nely, l1 = lmid; else l2 = lmid; end | Choose the next appropriate bounds for the bisection method |
| 81    | change = max(abs(xnew(:)-x(:))); | Calculate the L<sup>&infin;</sup> norm |
| 82    | x = xnew; | Prepare x for the next iteration |
| 84-85 | fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c, ...<br />&nbsp;&nbsp;mean(xPhys(:)),change); | Print the iteration number, objective function (compliance) value, volume fraction (this should be unchanged throughout the optimisation process), and change = L<sup>&infin;</sup> norm |
| 87    | colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow; | Plot the design (xPhys matrix) as a greyscale image |
| 88    | end | End function |
# Example of image insertion
![A test image](docs/images/example_drawing.svg)
*image_caption*

More text goes here ...
