function [poissFunc] = setup_poisson(pts,tri,epsCoeff,opt)
% Initializes the solver for the Poisson equation.
% 
% INPUTS:
%	pts: node locations, each column stores the coordinates for a particular node
%	tri: triangulation, each column stores the node-indices for a particular simplex
%	epsCoeff: coefficient representing the scalar/tensor field
%	method (optional): method to use for the poisson solver
%
% OUTPUTS:
%	poissFunc: returns the coefficients of the scalar field that solves the
%	Poisson equation. Use the syntax:
%		- v = poissFunc("rho",rhoCoeff):	for a known charge density
%		- v = q*poissFunc("psi",psiCoeff):	for a set of wavefunctions with charge q
%
% SEE ALSO: FEMAT, SETUP_LINSOLVE, SETUP_LITSOLVE
arguments
	pts (:,:) double {mustBeFinite,mustBeReal}
	tri (:,:) {mustBeInteger,mustBePositive}
	epsCoeff (:,:) double {mustBeFinite,mustBeReal}
	opt.method string {mustBeMember(opt.method,...
		["decompose","preinvert","full","iterative",...
		"lsqr","pcg","minres","gmres","qmr","tfqmr"])} = "decompose"
	opt.useZeroMean (1,1) logical = false;
end

% find the number of vertices per simplex
[Ndim,Npts] = size(pts);
tri = tri(1:Ndim+1,:);
[Nvtx,Ntri] = size(tri);
svol = simplex_vol(pts,tri);

% compute the left-hand side of the poisson equation
L = FEmat(pts,tri,"stiffness",epsCoeff);

% Poisson boundary matrix
if opt.useZeroMean == false
	sL = FEmat(pts,tri,"surface", epsCoeff);
	L = L - sL;
else
	% use zero-mean to enforce unique solution
	M = FEmat(pts,tri,"overlap");
	L = [[L, sum(M,2)];
		[sum(M,1), 0]];
end
L = L/(4*pi);

% select solver method
if ismember(opt.method,["decompose","preinvert","full"])
	solverFunc = setup_linsolve(L,1,opt.method);
elseif opt.method~="iterative"
	solverFunc = setup_litsolve(L,1,opt.method);
else
	solverFunc = setup_litsolve(L,1,"auto");
end

% rho is first-order on the domain (if rho has one argument)
[intCoeff,aa,bb] = bary_integral(Ndim,2);
int_ab = intCoeff .* svol;
t_aa = tri(aa,:);
t_bb = tri(bb,:);

% flatten for accumarray
int_ab = int_ab(:);
t_aa = t_aa(:);
t_bb = t_bb(:);

% rho is second-order on the domain (if rho has two arguments)
[intCoeff,ii,jj,kk] = bary_integral(Ndim,3);
int_ijk = intCoeff .* svol;
t_ii = tri(ii,:);
t_jj = tri(jj,:);
t_kk = tri(kk,:);

% flatten for accumarray
int_ijk = int_ijk(:);
t_ii = t_ii(:);
t_jj = t_jj(:);
t_kk = t_kk(:);

% return the density-vector constructor function
poissFunc = @compute_potential;

%-----------------------------------------------------------
	function [out] = rvec_from_rho(rho)
	% computes the right side of the Poisson equation from the density rho
	if isvector(rho) && numel(rho)==Npts
		rho = rho(t_bb);
	elseif isequal(size(rho),[Nvtx,Ntri])
		rho = rho(bb,:);
	elseif isequal(size(rho),[Ntri,Nvtx])
		rho = rho(:,bb).';
	else
		error("ERROR: the input rho has invalid size")
	end
	% integrate
	out = accumarray(t_aa, (int_ab .* rho(:)), [Npts,1]);
	end
%-----------------------------------------------------------
	function [out] = rvec_from_psi(psi)
	% computes the right side of the Poisson equation from the wavefunctions psi
	if ~ismatrix(psi)
		error("ERROR: the input psi has invalid size")
	elseif size(psi,1)==Npts
		rho = zeros(size(t_ii));
		for nn = 1:size(psi,2)
			ket = psi(:,nn);
			rho = rho + (conj(ket(t_jj)) .* ket(t_kk));
		end
	elseif size(psi,2)==Npts
		rho = zeros(size(t_ii));
		for nn = 1:size(psi,1)
			ket = psi(nn,:);
			rho = rho + (conj(ket(t_jj)) .* ket(t_kk));
		end
	elseif isequal(size(psi),[Nvtx,Ntri])
		rho = (conj(psi(jj,:)) .* psi(kk,:));
	elseif isequal(size(psi),[Ntri,Nvtx])
		rho = (conj(psi(:,jj)) .* psi(:,kk)).';
	else
		error("ERROR: the input psi has invalid size")
	end
	% integrate
	out = accumarray(t_ii, (int_ijk .* rho(:)), [Npts,1]);
	end
%-----------------------------------------------------------
	function [v_ee] = compute_potential(data_type, data_val)
	% returns the Coulombic potential in vector form
	if ~ismember(data_type,["rho","psi"])
		error("ERROR: the first argument must be either 'rho' or 'psi'.")
	elseif data_type=="rho"
		r = rvec_from_rho(data_val);
	elseif data_type=="psi"
		r = rvec_from_psi(data_val);
	end
	% switch between natural and zero-mean boundary conditions
	if opt.useZeroMean==false
		v_ee = solverFunc(r);
	else
		r = [r; 0];
		v_ee = solverFunc(r);
		v_ee = v_ee(1:end-1);
	end
	%
	end
%-----------------------------------------------------------
end