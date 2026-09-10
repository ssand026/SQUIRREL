function [pts,tri,edg] = meshify_geom(geom, numNodes, opt)
% Generate a mesh from a geometry description.
% The resulting mesh will approximately contain the requested number of interior nodes.
%
% REQUIRES:
%	Partial Differential Equation Toolbox
%
% SEE ALSO: DECSG, INITMESH
arguments
	geom
	numNodes (1,1) double {mustBeInteger,mustBePositive}
	opt.tol (1,1) double {mustBeInteger,mustBeInRange(opt.tol,0,numNodes)} = 2; 
	opt.maxIter (1,1) double {mustBeInteger,mustBePositive} = 128;
end

% get the total area of the geometry
[pts,~,tri] = initmesh(geom,"Hmax",Inf,"Jiggle","off","MesherVersion",'R2013a');
geomArea = sum(pdetrg(pts,tri));

% total area per node (with densest packing)
node_area = pi/(2*sqrt(3)) * geomArea/numNodes;

% get the average distance between nodes
node_dist = 2*sqrt(node_area/pi);
hmax = 1.2*node_dist;

z = 1.01; % scale factor for updating guesses
last = numNodes;
for iter = 1:opt.maxIter
	[pts,edg,tri] = initmesh(geom,"Hmax",hmax,"Jiggle","on","MesherVersion","R2013a");
	ext_nodes = unique(edg([1 2], any(edg([6 7],:)==0,1)));
	nodes = size(pts,2) - length(ext_nodes);
	
	if nodes < (numNodes - opt.tol)
		% hmax too big
		if last > numNodes; z = sqrt(z); end
		hmax = hmax/z;
	elseif nodes > (numNodes + opt.tol)
		% hmax too smol
		if last < numNodes; z = sqrt(z); end
		hmax = hmax*z;
	else
		break;
	end
	last = nodes;
	if z-1 <= 1e-5; break; end
end
end