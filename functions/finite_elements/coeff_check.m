function [coeff,order,continuity] = coeff_check(pts,tri,coeff,opt)
% Checks the dimensions of the inputted scalar/vector/tensor field
%
% INPUTS:
%	pts: node locations, each column stores the coordinates for a particular node
%	tri: triangulation, each column stores the node-indices for a particular simplex
%	coeff: numeric array representing a scalar/vector/tensor field
%	strictDims(optional,bool): require coefficient to be a pure field (no extra dims)
%
% OUTPUTS:
%	1) coeff: the re-sized input array
%	2) order: the field order i.e., scalar/vector/tensor 
%		- 0:	coeff is a scalar field, [d1,d2]==[1,1]
%		- 1:	coeff is a vector field, [d1,d2]==[Ndim,1]
%		- 2:	coeff is a tensor field, [d1,d2]==[Ndim,Ndim]
%	3) continuity: whether coeff is constant/continuous over the domain/simplexes
%		- "domain_constant":	size(coeff) = [d1,d2]
%		- "domain_continuous":	size(coeff) = [d1,d2,Npts]
%		- "simplex_constant":	size(coeff) = [d1,d2,Ntri]
%		- "simplex_continuous":	size(coeff) = [d1,d2,Nvtx,Ntri]
%
% SEE ALSO: COEFF_EVAL, COEFF_EXPAND, FUNC_REDUCE
arguments
	pts (:,:) double {mustBeFinite,mustBeReal}
	tri (:,:) {mustBeInteger,mustBePositive}
	coeff double {mustBeFinite,mustBeReal}
	opt.strictDims (1,1) logical = true;
end

% trim the triangulation to the nearest polynomial order
Ndim = size(pts,1);
vtx_len = cumsum(factorial(Ndim+1)./(factorial(Ndim-[0:Ndim]) .* factorial([1:Ndim+1])));
Npoly = find(size(tri,1)>=vtx_len,1,"last");
tri = tri(1:vtx_len(Npoly),:);

% check size of the coefficients
Ndim = size(pts,1); % number of spatial dimensions
Npts = size(pts,2); % number of nodes in the mesh
Nvtx = size(tri,1); % number of vertices per simplex
Ntri = size(tri,2); % number of simplices in the mesh

csize = size(coeff,[1:ndims(coeff)]);
cdims = (csize == [Ndim; Npts; Nvtx; Ntri]);

if any(sum(cdims,1) > 1)
	% dimensions coincide, cannot determine coefficient shape
	error("ERROR: overlap between the number of nodes, dimensions, vertices, or simplexes " + ...
		"of the mesh. These numbers are typically distinct outside of extremely small meshes")
end

% find which dimensions of coeff correspond to the various mesh measures (Ndim, Npts, Nvtx, Ntri)
dim_is_Ndim = find(cdims(1,:));
dim_is_Npts = find(cdims(2,:));
dim_is_Nvtx = find(cdims(3,:));
dim_is_Ntri = find(cdims(4,:));

% find any dimensions that do not correspond to the mesh measures
other_dim = find(~any(cdims,1));
other_dim = [other_dim(csize(other_dim)>1), other_dim(csize(other_dim)<=1)];

if opt.strictDims && any(csize(other_dim)~=1)
	error("ERROR: one or more of the array dimensions have a different size than the number of " + ...
		"nodes, dimensions, vertices, or simplexes in the mesh.")
end

% check the order/degree of the tensor field
% order = 0; scalar field
% order = 1; vector field
% order = 2; tensor field (matrices)
% order > 2; higher-order tensor field
order = numel(dim_is_Ndim);
if order > 2
	warning("Tensor coefficients above 2nd order are not currently supported")
end

% add check for whether the coefficient can be reduced to a more compact form ?

% check the coefficient continuity
dim_counts = [numel(dim_is_Npts), numel(dim_is_Nvtx), numel(dim_is_Ntri)];
if isequal(dim_counts,[0 0 0])
	% field is constant on the domain
	continuity = "domain_constant";
elseif isequal(dim_counts,[1 0 0])
	% field is continuous on the domain
	continuity = "domain_continuous";
elseif isequal(dim_counts,[0 1 0])
	% field is constant over each simplex
	continuity = "simplex_constant";
elseif isequal(dim_counts,[0 1 1])
	% field is continuous over each simplex
	continuity = "simplex_continuous";
else
	% field has dimensions incompatible with the allowed forms
	error("ERROR: the inputted coefficient has invalid shape")
end

% permute the coefficient dimensions
dimOrder = [dim_is_Ndim, dim_is_Npts, dim_is_Nvtx, dim_is_Ntri, other_dim];
dimOrder = [dimOrder, setdiff([1:ndims(coeff)],dimOrder)];
coeff = permute(coeff,dimOrder);
end