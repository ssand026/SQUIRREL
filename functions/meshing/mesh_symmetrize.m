function [pts,tri,edg,varargout] = mesh_symmetrize(pts,tri,edg,num_symm)
% Construct a fully symmetric mesh from a smaller sub-section.
% The input num_symm indicates the number of sub-regions required to create the
% full mesh, with num_symm==2 resulting in mirror symmetry across the x-axis.
arguments
	pts (:,:) double {mustBeFinite,mustBeReal}
	tri (:,:) {mustBeInteger,mustBePositive}
	edg (:,:) {mustBeNumeric,mustBeFinite}
	num_symm (1,1) double {mustBeInteger,mustBePositive}
end

[Ndim,Npts] = size(pts);
[pt_conn,~,~,~] = mesh_conn(pts,tri);

% keep track of sub-region indices
tri = [tri; ones(1,size(tri,2))];

% find edges that lie on the axes of symmetry
%-----------------------------------------------------------
do_reflect = (mod(num_symm,2)==0);
if do_reflect
	l_edge = at_angle(pts,edg,270);
	r_edge = at_angle(pts,edg,270+360/num_symm);
else
	l_edge = at_angle(pts,edg,270-360/num_symm);
	r_edge = at_angle(pts,edg,270+360/num_symm);
end

% convert to logical indicies
l_edge = any(edg(5,:)==l_edge,1);
r_edge = any(edg(5,:)==r_edge,1);
r_edge = r_edge & ~l_edge;

% find optimal sort order for the mesh
l_pts = unique(reshape(edg([1 2],l_edge,1),[],1));
r_pts = unique(reshape(edg([1 2],r_edge,1),[],1));
midpt = reshape(setdiff((1:Npts),[l_pts; r_pts]),[],1);

ord_l = l_pts(symrcm(pt_conn(l_pts,l_pts)));
ord_r = r_pts(symrcm(pt_conn(r_pts,r_pts)));
ord_m = midpt(symrcm(pt_conn(midpt,midpt)));

% delete redundant edges
edg(:, (l_edge | r_edge)) = [];

% re-order mesh
[pts,tri,edg] = mesh_reorder([ord_l; ord_m; ord_r],pts,tri,edg);

% reflect sub-mesh across the x-axis
%-----------------------------------------------------------
if do_reflect
	% get edges continous across x reflection
	x1 = pts(1,edg(1,:)); x2 = pts(1,edg(2,:));
	y1 = pts(2,edg(1,:)); y2 = pts(2,edg(2,:));
	flat_edge = [edg(5,:); y1==y2; x1==0 | x2==0];
	flat_edge = unique(flat_edge.','rows');
	flat_edge = flat_edge(flat_edge(:,2) & flat_edge(:,3),1);
	flat_indx = any(edg(5,:)==flat_edge,1);

	% make sure any flat edges are continuous upon reflection
	ie = @(list) full(sparse(list,1,1,7,1));
	f = flat_indx/2;
	edg(5,:) = num_symm * edg(5,:);
	e_l = edg .* ((ie([3 4]).*(~f + f)) + ~ie([3 4])) + (ie([3 4]) .* f);
	e_r = edg .* ((ie([3 4]).*(~f - f)) + ~ie([3 4])) + (ie([3 4]) .* f);

	% reflect across x=0
	Npts = size(pts,2);
	pts = [pts, [-1; 1].*pts];
	edg = [e_l, e_r - (ie(5) .* ~flat_indx) + (Npts * ie([1 2]))];
	tri = [tri, tri+[Npts; Npts; Npts; 0; 1]];
end

% rotate the sub-mesh to generate the full symmetry
%-----------------------------------------------------------
nssr = 1 + do_reflect; % total number of sub-symmetric regions
num_symm = num_symm / nssr;          % S indicates the number of remaining rotations
if num_symm > 1
	[Npts,ne,nt] = deal(size(pts,2),size(edg,2),size(tri,2));
	[pci,eci,tci] = deal([1:Npts],[1:ne],[1:nt]);
	pts = [pts, zeros([1 num_symm-1].*size(pts))];
	edg = [edg, zeros([1 num_symm-1].*size(edg))];
	tri = [tri, zeros([1 num_symm-1].*size(tri))];
	for kk = 1:(num_symm-1)
		pts(:,kk*Npts+pci) = exact_rotate(pts(:,pci), (-2 * kk / num_symm));
		edg(:,kk*ne+eci) = edg(:,eci) + kk * [Npts; Npts;  0; 0; -1; 0; 0];
		tri(:,kk*nt+tci) = tri(:,tci) + kk * [Npts; Npts; Npts; 0; nssr];
	end
end

% find shortest edge length
%-----------------------------------------------------------
[edg1,edg2] = ncombsk(Ndim+1,2);
min_edge = Inf;
for kk = 1:numel(edg1 & edg2)
	edge_lengths = vecnorm(pts(:,tri(edg1(kk),:)) - pts(:,tri(edg2(kk),:)));
	min_edge = min(min_edge,min(edge_lengths));
end

% find overlapping nodes and remove duplicates
%-----------------------------------------------------------
% find overlapping nodes
ptol = round(2*pts/min_edge).';
[~,ia,ib] = unique(ptol,"rows","stable");
bb = (1:length(ib)).';

% remove duplicate nodes
pts = pts(:,ia);

% remove duplicate edges
for kk = 1:Ndim
	new_ind = sum((edg(kk,:)==bb).*ib,1);
	edg(kk,any(edg(kk,:)==bb,1)) = 0;
	edg(kk,:) = edg(kk,:) + new_ind;
end
edg(:,edg(6,:)==edg(7,:))=[];

% remove duplicate triangles
for kk = 1:Ndim+1
	new_ind = sum((tri(kk,:)==bb).*ib,1);
	tri(kk,any(tri(kk,:)==bb,1)) = 0;
	tri(kk,:) = tri(kk,:) + new_ind;
end
tri(:, (tri(1,:)==tri(2,:) | tri(2,:)==tri(3,:) | tri(3,:)==tri(1,:)) ) = [];

% return the symmetry index
varargout{1} = tri(end,:);
tri = [sort(tri(1:Ndim+1,:),1); tri(Ndim+2,:)];
tri = tri(1:Ndim+2,:);

% re-index edges appropriately
edges = unique(edg(5,:));
edg = edg .* [1 1 1 1 1i 1 1].';
for kk = 1:length(edges)
	edg(5,edg(5,:)==1i*edges(kk)) = kk;
end
edg = real(edg);
edg = sortrows([edg([5 3 4],:).',edg.']).';  edg([1 2 3],:) = [];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [edges] = at_angle(p,e,theta)
% find edges on the symmetry axis specified by theta
x = p(1,:); y = p(2,:);
angle = [mod(atan2d(y(e(1,:))-y(e(2,:)),x(e(1,:))-x(e(2,:))),180);
	     mod(atan2d(y(e(2,:))+y(e(1,:)),x(e(2,:))+x(e(1,:))),360)];
tol = 10*max(eps(angle),[],2);
window = [mod(theta + tol(1)*[-1, 1] ,180);
	      mod(theta + tol(2)*[-1, 1] ,360)];
edges = [e(5,:); all(angle>=window(:,1) & angle<=window(:,2),1)];
edges = unique(edges.','rows');
edges = edges(~~edges(:,2),1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [out] = exact_rotate(p,theta)
% rotates points by given angle while maintaining exact precision
x0 = p(1,:);
y0 = p(2,:);
out(1,:) = cospi(theta) * x0 - sinpi(theta) * y0;
out(2,:) = sinpi(theta) * x0 + cospi(theta) * y0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
