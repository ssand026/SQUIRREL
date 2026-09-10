function [pts,tri,edg] = mesh_reorder(order,pts,tri,edg)
% Re-orders the mesh according to the specified node ordering
arguments
	order (:,1) {mustBeInteger,mustBePositive}
	pts (:,:) double {mustBeFinite,mustBeReal}
	tri (:,:) {mustBeInteger,mustBePositive}
	edg (:,:) {mustBeReal} = [];
end
[Ndim,Npts] = size(pts);

if any(order > Npts)
	error("ERROR: the specified ordering has indices outside the size of the mesh")
elseif numel(setdiff([1:Npts],order))>=1
	error("ERROR: specified ordering is not a complete re-arrangement of the original data")
end

% reorder vals in pts
pts = pts(:,order);

% list where each element contains the corresponding node's new location
[~,orderIndex] = sort(order(:));

% swap indices in tri
tri = [orderIndex(tri(1:Ndim+1,:)); tri(Ndim+2:end,:)];

% swap indices in edg
if ~isempty(edg)
	edg = [orderIndex(edg(1:Ndim,:)); edg(Ndim+1:end,:)];
end
end