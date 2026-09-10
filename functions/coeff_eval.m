function [c_vals,c_dims,c_type] = coeff_eval(pts,tri,coeff,reg)
% Evaluates the input as a finite-element coefficient.
%
% Takes a function_handle, cell, string, or numerical array and attempts to
% evaluate it as a scalar/vector/tensor field over the domain
arguments
	pts (:,:) double {mustBeFinite,mustBeReal}
	tri (:,:) {mustBeInteger,mustBePositive}
	coeff {mustBeA(coeff,{'function_handle','char','string','numeric','logical','cell'})}
	reg (:,:) {mustBeA(reg,{'logical','numeric'})} = []
end
% local type-sorting functions
isfunc = @(a) isa(a,"char") | isa(a,"string") | isa(a,"function_handle");
isnumb = @(a) isnumeric(a) | islogical(a);

% find the number of vertices per simplex
Ndim = size(pts,1);
vtx_len = cumsum(factorial(Ndim+1)./(factorial(Ndim-[0:Ndim]) .* factorial([1:Ndim+1])));
Npoly = find(size(tri,1)>=vtx_len,1,"last");

% trim triangulation if necessary
Nvtx = vtx_len(Npoly);
if size(tri,1)==Nvtx+1 & isempty(reg)
	reg = tri(Nvtx+1,:);
end
tri = tri(1:Nvtx,:);

% define the pre-set variables for position, mesh, etc.
norm_pts = pts./(max(max(pts,[],2)-min(pts,[],2))/2);

def_var = {};
if Ndim >= 1; def_var{end+1} = {"x",pts(1,:)}; def_var{end+1} = {"xx",norm_pts(1,:)}; end
if Ndim >= 2; def_var{end+1} = {"y",pts(2,:)}; def_var{end+1} = {"yy",norm_pts(2,:)}; end
if Ndim >= 3; def_var{end+1} = {"z",pts(3,:)}; def_var{end+1} = {"zz",norm_pts(3,:)}; end
if Ndim >= 2; def_var{end+1} = {"r",vecnorm(pts)}; def_var{end+1} = {"rr",vecnorm(norm_pts)}; end

if ~isempty(pts); def_var{end+1} = {"pts",pts}; end
if ~isempty(tri); def_var{end+1} = {"tri",tri}; end
if ~isempty(reg); def_var{end+1} = {"reg",reg}; end

% evaluate the coefficient
if isa(coeff,"cell")
	% cell: try to evaluate each cell as an individual coefficient
	c_vals = cell(size(coeff));
	c_dims = zeros(size(coeff));
	c_type = strings(size(coeff));
	
	for ii = 1:numel(coeff)
		if isfunc(coeff{ii})
			% string, char, or function
			[val,~,~,bool] = func_reduce(coeff{ii},def_var);
			if bool
				[c_vals{ii},c_dims(ii),c_type(ii)] = coeff_check(pts,tri,val);
			else
				c_vals{ii} = NaN;
				c_type(ii) = "no eval";
			end
		elseif isnumb(coeff{ii})
			% numeric
			[c_vals{ii},c_dims(ii),c_type(ii)] = coeff_check(pts,tri,coeff{ii});
		else
			c_type(ii) = "invalid type";
		end
	end
elseif isfunc(coeff)
	% string, char, or function: try evaluating the input, then check coefficient shape
	[val,~,~,bool] = func_reduce(coeff,def_var);
	if bool
		[c_vals,c_dims,c_type] = coeff_check(pts,tri,val);
	else
		c_vals = NaN;
		c_type = "no eval";
	end
elseif isnumb(coeff)
	% numeric type: check coefficient shape
	[c_vals,c_dims,c_type] = coeff_check(pts,tri,coeff);
else
	c_type = "invalid type";
end
end