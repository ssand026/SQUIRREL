function [pts,tri,edg] = mesh_symrcm(pts,tri,edg)
% Re-orders the mesh to minimize the bandwidth of the adjacency matrix.
% Will generally improve the performance of calculations with large/sparse systems.
%
% SEE ALSO: MESH_REORDER, SYMRCM
arguments
	pts (:,:) {mustBeFinite,mustBeReal}
	tri (:,:) {mustBeInteger,mustBePositive}
	edg (:,:) {mustBeFinite,mustBeReal} = [];
end
% get the new mesh order
adjM = mesh_adjacency(pts,tri);
order = symrcm(adjM);
[pts,tri,edg] = mesh_reorder(order,pts,tri,edg);
end