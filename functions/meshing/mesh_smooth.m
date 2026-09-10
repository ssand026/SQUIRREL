function [pts,tri] = mesh_smooth(pts,tri,iters)
% Improves the quality of a mesh by shifting each node towards the center of its connected neighbors
% The input "iters" determines how many times the node positions are re-adjusted,
% with larger values resulting in more smoothing of the mesh.
arguments
	pts (:,:) double {mustBeFinite,mustBeReal}
	tri (:,:) {mustBeInteger,mustBePositive}
	iters (1,1) {mustBeInteger,mustBePositive} = 1;
end

% logical array for nodes on the domain boundary
bnd_indx = boundary_nodes(pts,tri);
bnd_pts = pts(:,bnd_indx);

% get the normalized adjacency matrix
adjM = mesh_adjacency(pts,tri);
adjM = adjM ./ sum(adjM,1);

% begin averaging node positions
for jj = 1:iters
	% apply the normalized adjacency matrix
	pts = (pts * adjM);

	% keep the boundary points fixed
	pts(:,bnd_indx) = bnd_pts;
end
% done
end