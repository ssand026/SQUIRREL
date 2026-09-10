function [psi_array] = prop_PA(psi0, field, H0, mu, tau, opt)
% Propagates states using a scaled and squared Pade approximant.
arguments
	psi0  (:,:) double {mustBeFinite}
	field (:,:) double {mustBeFinite}
	H0    (:,:) double {validateattributes(H0,{'double'},{'finite','square'})}
	mu    (1,:) cell
	tau   (1,1) double {mustBeFinite}
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
		psi_array = prop_PA_fwd(psi_array, field, H0, mu, tau);
	case "bwd"
		% initial states are row vectors
		psi_array = complex(zeros(Nrow,Ncol,Ntau+1));
		psi_array(:,:,Ntau+1) = psi0;
		psi_array = prop_PA_bwd(psi_array, field, H0, mu, tau);
end
% done
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [psi_array] = prop_PA_fwd(psi_array, field, H0, mu, tau)
% forwards-propagate a set of states using Pade Approximation
[Ndim, Ntau] = size(field);
psi = psi_array(:,:,1);
orig_norm = sum((abs(psi).^2),1);

for jj = 1:Ntau
	% time-dependent Hamiltonian
	Ht = H0;
	for kk = 1:Ndim
		Ht = Ht + (field(kk,jj) * mu{kk});
	end
	% propagate
	psi = expm(-1i * tau * Ht) * psi;
	% normalize and save
	curr_norm = sum((abs(psi).^2),1);
	psi = psi .* sqrt(orig_norm ./ curr_norm);
	psi_array(:,:,jj+1) = psi;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [psi_array] = prop_PA_bwd(psi_array, field, H0, mu, tau)
% backwards-propagate a set of states using Pade Approximation
[Ndim, Ntau] = size(field);
psi = psi_array(:,:,Ntau+1);
orig_norm = sum((abs(psi).^2),2);

for jj = flip(1:Ntau)
	% time-dependent Hamiltonian
	Ht = H0;
	for kk = 1:Ndim
		Ht = Ht + (field(kk,jj) * mu{kk});
	end
	% propagate
	psi = psi * expm(-1i * tau * Ht);
	% normalize and save
	curr_norm = sum((abs(psi).^2),2);
	psi = psi .* sqrt(orig_norm ./ curr_norm);
	psi_array(:,:,jj) = psi;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%