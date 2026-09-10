function [pts,tri,edg] = mesh_coreshell(side_lens, num_sides, num_nodes, symmetry)
% generates a mesh for a 2D core-shell nanowire
arguments
	side_lens (1,:) double {mustBePositive}
	num_sides (1,1) double {mustBeGreaterThan(num_sides,2)}
	num_nodes (1,1) double {mustBePositive}
	symmetry string {mustBeMember(symmetry,["none","half","full"])}
end

[geom,s] = coreshell_geom(side_lens, num_sides, symmetry);
[pts,tri,edg] = meshify_geom(geom,ceil(num_nodes/s));
[pts,tri] = mesh_smooth(pts,tri,10);
[pts,tri] = mesh_symrcm(pts,tri,edg);

[pts,tri,edg,sindx] = mesh_symmetrize(pts,tri,edg,s);
[pts,tri] = mesh_anneal(pts,tri,"symmIndex",sindx);
[pts,tri] = mesh_smooth(pts,tri,10);

valid_edgs = (edg(6,:)~=edg(7,:)) & (edg(1,:)~=edg(2,:));
edg = edg(:,valid_edgs);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [geom,s] = coreshell_geom(side_lengths, num_sides, symmetry)
% generates the geometry matrix for a polygonal core-shell geometry
side_lengths = sort(side_lengths(:));
num_layers = numel(side_lengths);
sub_pts = cell(1,num_layers);

switch symmetry
	case "none"
		s = 1;
		for jj = 1:num_layers
			sub_pts{jj} = polygon_full(num_sides,side_lengths(jj));
		end
		geom = full_geom(sub_pts);
	case "half"
		s = 2;
		for jj = 1:num_layers
			sub_pts{jj} = polygon_half(num_sides,side_lengths(jj));
		end
		geom = slice_geom(sub_pts);
	case "full"
		s = 2*num_sides;
		for jj = 1:num_layers
			sub_pts{jj} = polygon_segment(num_sides,side_lengths(jj));
		end
		geom = slice_geom(sub_pts);
	otherwise
		error("ERROR: invalid symmetry specification")
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [out] = polygon_segment(num_sides,side_len)
% minimum symmetrical segment of a polygonal region
r = side_len/(2 * sinpi(1/num_sides));
theta = 3/2 - 1/num_sides;
out = r * [0 cospi(theta); sinpi(theta) sinpi(theta)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [out] = polygon_half(num_sides,side_len)
% minimum symmetrical segment of a polygonal region
r = side_len/(2 * sinpi(1/num_sides));
theta = 3/2 - 1/num_sides;

len = floor(num_sides/2);
x = zeros(1,len+2);
y = zeros(1,len+2);
for jj = 0:len+1
	x(jj+1) = cospi(theta + 2*jj/num_sides);
	y(jj+1) = sinpi(theta + 2*jj/num_sides);
end
x(1)=0; x(end) = 0;
out = r*[x; y];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [out] = polygon_full(num_sides,side_len)
% full set of points for a polygonal region
r = side_len/(2 * sinpi(1/num_sides));
theta = 3/2 - 1/num_sides;

x = zeros(1,num_sides+1);
y = zeros(1,num_sides+1);
for jj = 0:num_sides
	x(jj+1) = cospi(theta + 2*jj/num_sides);
	y(jj+1) = sinpi(theta + 2*jj/num_sides);
end
out = r*[x;y];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [geom] = slice_geom(shapes)
% assembles geometry description matrix from a slice
num_shapes = max(size(shapes));
[xi,xf,yi,yf] = deal(0,0,0,0);

geom = [];
for k = 1:num_shapes
    pos = shapes{k};
    num_pts = size(pos,2)+1;
    
    g = ones(7,num_pts);
    g(1,:) = 2;
    g(2,:) = [xi, pos(1,:)];
    g(3,:) = [pos(1,:), xf];
    g(4,:) = [yi, pos(2,:)];
    g(5,:) = [pos(2,:), yf];
    g(6,:) = k;
    g(7,:) = k+1;
	g(7,[1 end]) = 0;
    if k==num_shapes; g(7,:) = 0; end
	[xi,xf,yi,yf] = deal(pos(1,1),pos(1,end),pos(2,1),pos(2,end));
	geom = [geom, g];
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [geom] = full_geom(shapes)
% assembles geometry description matrix for a full polygon
num_shapes = max(size(shapes));

geom = [];
for k = 1:num_shapes
    pos = shapes{k};
    num_pts = size(pos,2)-1;
    
    g = ones(7,num_pts);
    g(1,:) = 2;
    g(2,:) = pos(1,1:(end-1));
    g(3,:) = pos(1,2:(end-0));
    g(4,:) = pos(2,1:(end-1));
    g(5,:) = pos(2,2:(end-0));
    g(6,:) = k;
    g(7,:) = k+1;
    if k==num_shapes; g(7,:) = 0; end
    geom = [geom, g];
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%