function [B] = femat_nullspace(pts,tri)
% (memoized) Returns the finite-element nullspace matrix
arguments
	pts (:,:) double {mustBeFinite,mustBeReal}
	tri (:,:) {mustBeInteger,mustBePositive}
end

% enable memoization for repeated calls
persistent memFunc
if isempty(memFunc)
	memFunc = memoize(@femat_nullspace_main);
	memFunc.CacheSize = 1;
end
B = memFunc(pts,tri);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [B] = femat_nullspace_main(pts,tri)
Npts = size(pts,2);
bnodes = boundary_nodes(pts,tri);
B = speye(Npts,Npts);
B(:,bnodes) = [];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%