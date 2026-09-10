function [out] = coeff_expand(pts,tri,coeff)
% Expands the coefficient so its value is given for the vertices of each simplex
%
% The resulting 4-D array will have dimensions [d1, d2, Nvtx, Ntri], where 
% (d1==d2==1) for a scalar field, (d1==Ndim) & (d2==1) for a vector field, and 
% (d1==d2==Ndim) for a tensor field.
arguments
	pts (:,:) double {mustBeFinite,mustBeReal}
	tri (:,:) {mustBeInteger,mustBePositive}
	coeff (:,:,:,:) double
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

% squeeze dimensions of length 1 from the coefficient
oneDims = (size(coeff,[1:ndims(coeff)])==1);
coeff = permute(coeff,[find(~oneDims),find(oneDims)]);

% find which dimensions of coeff correspond to the various mesh measures (Ndim, Npts, Nvtx, Ntri)
csize = size(coeff,[1:ndims(coeff)]);
dim_is_Ndim = find(csize==Ndim & csize~=1);
dim_is_Npts = find(csize==Npts);
dim_is_Nvtx = find(csize==Nvtx);
dim_is_Ntri = find(csize==Ntri);

% re-order the coefficient dimensions
allDims = [1:ndims(coeff)];
dimOrder = [dim_is_Ndim, dim_is_Npts, dim_is_Nvtx, dim_is_Ntri];
dimOrder = [dimOrder, allDims(~ismember(allDims,dimOrder))];
coeff = permute(coeff,dimOrder);

% get the tensor-order of the coefficient
order = numel(dim_is_Ndim);
if order > 2
	error("ERROR: the coefficient cannot be a tensor of order 3 or higher.")
end

% reshape the coefficient into a vertex-based tensor
dim_counts = [numel(dim_is_Npts), numel(dim_is_Nvtx), numel(dim_is_Ntri)];

if isequal(dim_counts,[0 0 0])
	% field is constant over the domain
	if order==0
		out = repmat(coeff,Nvtx,Ntri);
	elseif order==1
		out = repmat(coeff,1,Nvtx,Ntri);
	elseif order==2
		out = repmat(coeff,1,1,Nvtx,Ntri);
	end
elseif isequal(dim_counts,[1 0 0])
	% field is continuous over the domain
	if order==0
		out = coeff(tri);
	elseif order==1
		out = zeros(Ndim,Nvtx,Ntri);
		for ii = 1:Ndim
			temp = reshape(coeff(ii,:),Npts,1);
			out(ii,:,:) = temp(tri);
		end
	elseif order==2
		out = zeros(Ndim,Ndim,Nvtx,Ntri);
		for ii = 1:Ndim
			for jj = 1:Ndim
				temp = reshape(coeff(ii,jj,:),Npts,1);
				out(ii,jj,:,:) = temp(tri);
			end
		end
	end
elseif isequal(dim_counts,[0 0 1])
	% field is constant over each simplex
	if order==0
		out = repmat(coeff.',Nvtx,1);
	elseif order==1
		out = repmat(coeff,1,Nvtx,1);
	elseif order==2
		out = repmat(coeff,1,1,Nvtx,1);
	end
elseif isequal(dim_counts,[0 1 1])
	% field is continuous over each simplex
	out = coeff;
else
	% field has dimensions incompatible with the allowed forms
	error("ERROR: the coefficient has invalid shape")
end

% final reshape step
if order==0
	out = permute(out,[3,4,1,2]);
elseif order==1
	out = permute(out,[1,4,2,3]);
elseif order==2
	% already in proper order
end
end