function [pts,tri,edgeLen] = meshify_STL(stl_file, numNodes)
% Generates a mesh from an STL file.
% The resulting mesh will approximately contain the requested number of nodes.
%
% REQUIRES:
%	Partial Differential Equation Toolbox
%
% SEE ALSO: STLREAD
arguments
	stl_file (1,1) string
	numNodes (1,1) double
end
% load geometry
geom = importGeometry(stl_file);
TR = stlread(stl_file);

% compute radius for close packing
[~,vol] = convhulln(TR.Points);
Ndim = size(TR.Points,2);
if Ndim==3
	close_pack = pi/(3*sqrt(2));
	svol = (vol * 0.95 * close_pack)/numNodes;
	r = (svol/(4/3*pi))^(1/3);
elseif Ndim==2
	close_pack = pi/(2*sqrt(3));
	svol = (vol * 0.95  * close_pack)/numNodes;
	r = (svol/pi)^(1/2);
end
Hgrad = sqrt(2);

% create mesh
[stlGeom] = makeMesh(geom,2*r,r,Hgrad);
pts = stlGeom.Nodes;
tri = stlGeom.Elements;

% normalize and center geometry
minima = min(pts,[],2);
maxima = max(pts,[],2);
scl = max(abs(maxima-minima));
pts = pts - (minima+maxima)/2;
pts = pts/scl;
edgeLen = 2*r / scl;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [geom_mesh] = makeMesh(geom,Hmax,Hmin,Hgrad)
% creates mesh from discreteGeometry object
arguments
	geom
	Hmax (1,1) double {mustBeNonnegative}
	Hmin (1,1) double {mustBeNonnegative}
	Hgrad (1,1) double {mustBeInRange(Hgrad,1,2)} = 4/3;
end
geomOrder = 'linear';

% sanity checks
if Hmax ~= 0 && Hmin ~= 0 && Hmin > Hmax
	error(message('pde:pdeModel:invalidHmaxHmin'));
end

% strip and check local refinement data
Hface = getValidLocalSizes([],geom.NumFaces, ...
	'pde:pdeGeometryClass:InvalidFaceIndex','pde:pdeModel:invalidHface');
Hedge = getValidLocalSizes([],geom.NumEdges, ...
	'pde:pdeGeometryClass:InvalidEdgeIndex','pde:pdeModel:invalidHedge');
Hvert = getValidLocalSizes([],geom.NumVertices, ...
	'pde:pdeGeometryClass:InvalidVertexIndex','pde:pdeModel:invalidHvertex');

% generate mesh
try
	[nodes,tet,cas,fas,eas,vas,Hmax,Hmin,Hgrad] = genmeshinternal(geom, ...
		Hmax,Hmin,Hgrad,Hface,Hedge,Hvert,geomOrder);
catch ex
	throwAsCaller(ex);
end
assoc = pde.FEMeshAssociation(tet, cas, fas, eas, vas);
geom_mesh = pde.FEMesh(nodes, tet, Hmax, Hmin, Hgrad, geomOrder, assoc);

%-----------------------------------------------------------
	function sizes = getValidLocalSizes(hc,num,err_index,err_format)
	sizes = zeros(1,num);
	if isempty(hc)
		return
	end
	if ~isvector(hc)
		error(message(err_format));
	end
	for i = 1:2:length(hc)
		id = hc{i};
		if (~isvector(id) || ~isreal(id) || ischar(id) || issparse(id(:)) || ...
				any(~isfinite(id(:))) || any(id < 1) || any(id > num))
			error(message(err_index));
		end
		sz = hc{i+1};
		if ~isreal(sz) || ~isscalar(sz) || ischar(sz) || sz < 0 || issparse(sz) || ~isfinite(sz)
			error(message('pde:pdeModel:invalidHlocal'));
		end
		sizes(id) = sz;
	end
	end
%-----------------------------------------------------------
end