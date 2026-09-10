function [adjM] = mesh_adjacency(pts,tri)
% Constructs the node adjacency matrix for the given mesh
arguments
	pts (:,:) double {mustBeFinite,mustBeReal}
	tri (:,:) {mustBeInteger,mustBePositive}
end
[Ndim,Npts] = size(pts);

% estimate the number of nonzero elements
Nnnz = (Ndim^2+Ndim)*Npts;

% get the adjacency matrix
adjM = logical(spalloc(Npts,Npts,Nnnz));
for aa = 1:Ndim+1
	for bb = aa+1:Ndim+1
		adjM = adjM | sparse(tri(aa,:),tri(bb,:),true,Npts,Npts);
	end
end

adjM = (adjM.' | adjM);
end