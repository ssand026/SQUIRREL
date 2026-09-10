function [convertFunc] = setup_operator(pts,tri)
% Returns a function that converts the scalar coefficient into its operator form.
arguments
	pts (:,:) double {mustBeFinite,mustBeReal}
	tri (:,:) {mustBeInteger,mustBePositive}
end
[Ndim,Npts] = size(pts);
[Nvtx,Ntri] = size(tri);

% compute integral coefficients
svol = simplex_vol(pts,tri);
[intCoeff,ii,jj,kk] = bary_integral(Ndim,3);
int_ijk = intCoeff .* svol;
t_ii = tri(ii,:);
t_jj = tri(jj,:);
t_kk = tri(kk,:);

% pass reference to the conversion function
convertFunc = @conversion_func;
%-----------------------------------------------------------
	function [V] = conversion_func(coeff)
	% reshape the input coeff
	coeff = squeeze(coeff);
	if isequal(size(coeff),[Npts,1]) || isequal(size(coeff),[1,Npts])
		coeff = coeff(t_kk);
	elseif isequal(size(coeff),[Ndim+1,Ntri])
		coeff = coeff(kk,:);
	elseif isequal(size(coeff),[Ntri,Ndim+1])
		coeff = coeff(:,kk).';
	else
		error("ERROR: input coefficient has invalid size")
	end
	% convert to operator
	Vij = sparse(t_ii, t_jj, int_ijk .* coeff, Npts, Npts);
	V = (Vij + Vij')/2;
	end
%-----------------------------------------------------------
end