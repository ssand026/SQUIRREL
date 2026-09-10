function [pts,tri,edg] = mesh_symamd(pts,tri,edg)
% Re-orders the mesh to maximize the sparsity of the adjacency matrix's Cholesky decomposition
% Will generally improve the performance of linear solves in large/sparse systems.
%
% SEE ALSO: SYMAMD
arguments
	pts (:,:) {mustBeFinite,mustBeReal}
	tri (:,:) {mustBeInteger,mustBePositive}
	edg (:,:) {mustBeFinite,mustBeReal} = [];
end
% get the new mesh order
adjM = mesh_adjacency(pts,tri);
order = symamd(adjM);
[pts,tri,edg] = mesh_reorder(order,pts,tri,edg);
end