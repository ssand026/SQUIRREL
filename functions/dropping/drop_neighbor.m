function [H] = drop_neighbor(H,pts,tri,max_edge)
% Drops elements from the Hamiltonian using a connection-based criterion
%
% Returns a dropped form of the matrix H, where the element H_ij -> 0 if the 
% number of edges connecting the ith and jth node is greater than "max_edge"
%
% SEE ALSO: DROP_PERTURB, DROP_DISTANCE, SPARSIFY
arguments
	H (:,:) double {mustBeFinite}
	pts (:,:) double {mustBeFinite,mustBeReal}
	tri (:,:) {mustBeInteger,mustBePositive}
	max_edge (1,1) {mustBeInteger,mustBePositive}
end

% exponentiate adjacency matrix to get the mask
adjM = mesh_adjacency(pts,tri);
to_drop = ~(speye(size(adjM)) | (adjM^(max_edge)));

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