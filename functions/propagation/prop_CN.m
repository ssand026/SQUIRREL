function [psi_array] = prop_CN(psi0, field, H0, mu, tau, M, opt)
% Propagates states using the Crank-Nicolson propagator.
arguments
	psi0  (:,:) double {mustBeFinite}
	field (:,:) double {mustBeFinite}
	H0    (:,:) double {validateattributes(H0,{'double'},{'finite','square'})}
	mu    (1,:) cell
	tau   (1,1) double {mustBeFinite}
	M     (:,:) double {validateattributes(M,{'double'},{'finite','square'})} = [];
	opt.direction string {mustBeMember(opt.direction,["fwd","bwd"])} = "fwd";
end
[Npts, ~] = size(H0);
[Ndim, Ntau] = size(field);
[Nrow, Ncol] = size(psi0);

if length(mu)~=Ndim
	% test if the number of dipole moment matrices equals the number of field dimensions
	error('ERROR: the number of field dimensions does not match mu')
elseif (opt.direction=="fwd") & Nrow~=Npts
	% test if the initial state has valid dimensions
	error('ERROR: the number of rows in the initial state do not match the Hamiltonian')
elseif (opt.direction=="bwd") & Ncol~=Npts
	% test if the initial state has valid dimensions
	error('ERROR: the number of columns in the initial state do not match the Hamiltonian')
end

switch opt.direction
	case "fwd"
		% initial states are column vectors
		psi_array = complex(zeros(Nrow,Ncol,Ntau+1));
		psi_array(:,:,1) = psi0;
		psi_array = prop_CN_fwd(psi_array, field, H0, mu, tau, M);
	case "bwd"
		% initial states are row vectors
		psi_array = complex(zeros(Nrow,Ncol,Ntau+1));
		psi_array(:,:,Ntau+1) = psi0;
		psi_array = prop_CN_bwd(psi_array, field, H0, mu, tau, M);
end
% done
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [psi_array] = prop_CN_fwd(psi_array, field, H0, mu, tau, M)
% forwards-propagate a set of states using the Crank-Nicolson method
[Ndim, Ntau] = size(field);
psi = psi_array(:,:,1);

if isempty(M); M = eye(size(H0),"like",H0); end
orig_norm = abs(sum(conj(psi).*(M*psi),1));

for jj = 1:Ntau
	% time-dependent Hamiltonian
	Ht = H0;
	for kk = 1:Ndim
		Ht = Ht + field(kk,jj) * mu{kk};
	end
	% propagate
	psi = (M - (1i*tau/2) * Ht) * psi;
	psi = (M + (1i*tau/2) * Ht) \ psi;
	% normalize and save
	curr_norm = abs(sum(conj(psi).*(M*psi),1));
	psi = psi .* sqrt(orig_norm ./ curr_norm);
	psi_array(:,:,jj+1) = psi;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [psi_array] = prop_CN_bwd(psi_array, field, H0, mu, tau, M)
% backwards-propagate a set of states using the Crank-Nicolson method
[Ndim, Ntau] = size(field);
psi = psi_array(:,:,Ntau+1);

if isempty(M); M = eye(size(H0),"like",H0); end
orig_norm = abs(sum((psi*M).*conj(psi),2));

for jj = flip(1:Ntau)
	% time-dependent Hamiltonian
	Ht = H0;
	for kk = 1:Ndim
		Ht = Ht + (field(kk,jj) * mu{kk});
	end
	% propagate
	psi = psi / (M + (1i*tau/2) * Ht);
	psi = psi * (M - (1i*tau/2) * Ht);
	% normalize and save
	curr_norm = abs(sum((psi*M).*conj(psi),2));
	psi = psi .* sqrt(orig_norm ./ curr_norm);
	psi_array(:,:,jj) = psi;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%