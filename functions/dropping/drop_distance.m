function [H] = drop_distance(H,pts,tri,max_dist)
% Drops elements from the Hamiltonian using a distance-based criterion
%
% Returns a dropped form of the  matrix H, where the element H_ij -> 0 if the
% distance between the ith and jth node is greater than "max_dist"
%
% SEE ALSO: DROP_PERTURB, DROP_NEIGHBOR, SPARSIFY
arguments
	H (:,:) double {mustBeFinite}
	pts (:,:) double {mustBeFinite,mustBeReal}
	tri (:,:) {mustBeInteger,mustBePositive}
	max_dist (1,1) {mustBeFinite,mustBePositive}
end
[Ndim,Npts] = size(pts);

% get the square of the distance between points
dist = zeros(Npts,Npts);
for jj = 1:Ndim
	rjj = meshgrid(pts(jj,:));
	dist = dist + (rjj-rjj.').^2;
end
to_drop = logical(dist > max_dist^2);

% check the size of the mask
B = femat_nullspace(pts,tri);
if all(size(H)==min(size(B))) && all(size(to_drop)==max(size(B)))
	% convert the mask to the reduced basis
	to_drop = logical(to_basis("sparse",to_drop,"mat",B));
elseif ~isequal(size(H),size(to_drop))
	error("ERROR: the inputted matrix is not compatible with the mesh size")
end

% apply the mask
H(to_drop) = 0;
end