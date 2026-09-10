function [B] = femat_orthobasis(pts,tri)
% (memoized) Returns the matrix that converts between the full and orthonormal basis
% 
% The orthobasis matrix is equal to the product of the nullspace matrix and the 
% square root of the overlap matrix.
arguments
	pts (:,:) double {mustBeFinite,mustBeReal}
	tri (:,:) {mustBeInteger,mustBePositive}
end

% enable memoization for repeated calls
persistent memFunc
if isempty(memFunc)
	memFunc = memoize(@femat_orthobasis_main);
	memFunc.CacheSize = 1;
end
B = memFunc(pts,tri);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [B] = femat_orthobasis_main(pts,tri)
% dirichlet nullspace-matrix
D = femat_nullspace(pts,tri);

% reduced-space mass matrix
M = femat_overlap(pts,tri);
M = (D' * M * D);

% orthonormal-basis conversion matrix
B = inv(sqrtm(full(M)));
B = (B + B.')/2;
B = full(D) * B;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%