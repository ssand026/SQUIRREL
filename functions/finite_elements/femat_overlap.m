function [M] = femat_overlap(pts,tri)
% (memoized) Returns the finite-element overlap matrix
arguments
	pts (:,:) double {mustBeFinite,mustBeReal}
	tri (:,:) {mustBeInteger,mustBePositive}
end

% enable memoization for repeated calls
persistent memFunc
if isempty(memFunc)
	memFunc = memoize(@femat_overlap_main);
	memFunc.CacheSize = 1;
end
M = memFunc(pts,tri);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [M] = femat_overlap_main(pts,tri)
% get the simplex volumes
[Ndim,Npts] = size(pts);
svol = simplex_vol(pts,tri);

% assemble integration coefficients
[intCoeff,ii,jj] = bary_integral(Ndim,2);

M = sparse(tri(ii,:),tri(jj,:),svol.*intCoeff,Npts,Npts);
M = (M + M.')/2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%