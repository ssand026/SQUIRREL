function [psi_array] = prop_CN_multi(psi0, field, H0, mu, tau, M, V_psi, opt)
% Propagates a set of interacting states using the Crank-Nicolson method.
arguments
	psi0  (:,:) double {mustBeFinite}
	field (:,:) double {mustBeFinite}
	H0    (:,:) double {validateattributes(H0,{'double'},{'finite','square'})}
	mu    (1,:) cell
	tau   (1,1) double {mustBeFinite}
	M     (:,:) double {validateattributes(M,{'double'},{'finite','square'})}
	V_psi  function_handle
	opt.direction string {mustBeMember(opt.direction,["fwd","bwd"])} = "fwd";
	opt.maxIters (1,1) {mustBePositive,mustBeInteger} = 12;
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

if isempty(M); M = eye(size(H0),"like",H0); end
switch opt.direction
	case "fwd"
		% initial states are column vectors
		psi_array = complex(zeros(Nrow,Ncol,Ntau+1));
		psi_array(:,:,1) = psi0;
		psi_array = prop_CN_fwd(psi_array, field, H0, mu, tau, M, V_psi, opt.maxIters);
	case "bwd"
		% initial states are row vectors
		psi_array = complex(zeros(Nrow,Ncol,Ntau+1));
		psi_array(:,:,Ntau+1) = psi0;
		psi_array = prop_CN_bwd(psi_array, field, H0, mu, tau, M, V_psi, opt.maxIters);
end
% done
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [psi_array] = prop_CN_fwd(psi_array, field, H0, mu, tau, M, V_psi, maxIters)
% forwards-propagate a set of states using the Crank-Nicolson method
[Ndim, Ntau] = size(field);
psi_a = psi_array(:,:,1);
orig_norm = abs(sum(conj(psi_a).*(M*psi_a),1));

for jj = 1:Ntau
	% time-dependent Hamiltonian
	Ht = H0;
	for kk = 1:Ndim
		Ht = Ht + field(kk,jj) * mu{kk};
	end
	
	% self-consistent propagation
	psi_next = psi_a;
	for kk = 1:maxIters
		% update Hamiltonian
		H_psi = Ht + V_psi((psi_a+psi_next)/2);
		
		% propagate
		psi_b = (M - (1i*tau/2) * H_psi) * psi_a;
		psi_c = (M + (1i*tau/2) * H_psi) \ psi_b;
		
		% check for convergence
		psi_diff = (psi_c - psi_next);
		diff_mag = abs(sum(conj(psi_diff).*(M*psi_diff),1));
		if all(diff_mag <= 2*eps)
			break
		else
			psi_next = psi_c;
		end
	end
	psi_a = psi_c;

	% normalize and save
	curr_norm = abs(sum(conj(psi_a).*(M*psi_a),1));
	psi_a = psi_a .* sqrt(orig_norm ./ curr_norm);
	psi_array(:,:,jj+1) = psi_a;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [psi_array] = prop_CN_bwd(psi_array, field, H0, mu, tau, M, V_psi, maxIters)
% backwards-propagate a set of states using the Crank-Nicolson method
[Ndim, Ntau] = size(field);
psi_c = psi_array(:,:,Ntau+1);
orig_norm = abs(sum((psi_c*M).*conj(psi_c),2));

for jj = flip(1:Ntau)
	% time-dependent Hamiltonian
	Ht = H0;
	for kk = 1:Ndim
		Ht = Ht + (field(kk,jj) * mu{kk});
	end

	% self-consistent propagation
	psi_next = psi_c;
	for kk = 1:maxIters
		% update Hamiltonian
		H_psi = Ht + V_psi((psi_c.'+psi_next.')/2);

		% propagate
		psi_b = psi_c / (M + (1i*tau/2) * H_psi);
		psi_a = psi_b * (M - (1i*tau/2) * H_psi);
		
		% check for convergence
		psi_diff = (psi_a - psi_next);
		diff_mag = abs(sum((psi_diff*M).*conj(psi_diff),2));
		if all(diff_mag <= 2*eps)
			break
		else
			psi_next = psi_a;
		end
	end
	psi_c = psi_a;

	% normalize and save
	curr_norm = abs(sum((psi_c*M).*conj(psi_c),2));
	psi_c = psi_c .* sqrt(orig_norm ./ curr_norm);
	psi_array(:,:,jj) = psi_c;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%