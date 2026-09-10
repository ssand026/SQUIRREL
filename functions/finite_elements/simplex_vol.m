function [vol] = simplex_vol(pts,tri)
% (memoized) computes the hypervolume of each simplex in the mesh
arguments
	pts (:,:) double {mustBeFinite,mustBeReal}
	tri (:,:) {mustBeInteger,mustBePositive}
end

% enable memoization for repeated calls
persistent memFunc
if isempty(memFunc)
	memFunc = memoize(@simplex_vol_main);
	memFunc.CacheSize = 1;
end
vol = memFunc(pts,tri);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vol] = simplex_vol_main(pts,tri)
Ndim = size(pts,1);
Ntri = size(tri,2);

T = pts(:,tri(1:Ndim+1,:)); T(Ndim+1,:) = 1;
T = reshape(T,Ndim+1,Ndim+1,Ntri);

detT = zeros(1,Ntri);
for ss = 1:Ntri
	detT(ss) = abs(det(T(:,:,ss)));
end
vol = detT / factorial(Ndim);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%