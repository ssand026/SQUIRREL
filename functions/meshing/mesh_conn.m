function [pt_conn, num_tri, bnd_indx, ext_indx] = mesh_conn(pts,tri)
% Finds the connectivity and edge/boundary indices for a 2D mesh
%
% For the edge with nodes i,j, there are two triangle indices n & m such that
% tri(:,m) & tri(:,n) both include i & j. The connectivity matrix pt_conn is
% configured such that pt_conn(i,j)==n & pt_conn(j,i)==m
arguments
	pts (2,:) double {mustBeFinite,mustBeReal}
	tri (:,:) {mustBeInteger,mustBePositive}
end
Npts = size(pts,2);
Ntri = size(tri,2);

% construct the connectivity matrix
connM = zeros(Npts,Npts);
for jj = 1:Ntri
	a = tri(1,jj); b = tri(2,jj); c = tri(3,jj);
	if (connM(a,b)==0); connM(a,b) = jj; elseif (connM(b,a)==0); connM(b,a) = jj; end
	if (connM(b,c)==0); connM(b,c) = jj; elseif (connM(c,b)==0); connM(c,b) = jj; end
	if (connM(c,a)==0); connM(c,a) = jj; elseif (connM(a,c)==0); connM(a,c) = jj; end
end
pt_conn = sparse(connM);
is_conn = logical(pt_conn);

% find edges with only a single connected triangle (external edges)
ext_indx = (is_conn ~= is_conn.');

% move any external edges to upper part of the matrix
upper_ext = triu((pt_conn + pt_conn.') .* ext_indx);
pt_conn(ext_indx) = 0;
pt_conn = pt_conn + upper_ext;

% edges on domain boundaries
[row,col,triIndx] = find(pt_conn);
bnd_indx = sparse(row,col,tri(4,triIndx),Npts,Npts);
bnd_indx = (bnd_indx ~= bnd_indx.');
bnd_indx = bnd_indx | ext_indx;

% number of connected triangles
num_tri = full(sum(is_conn,2)-any(ext_indx,2));
end