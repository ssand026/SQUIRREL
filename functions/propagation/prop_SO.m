function [psi_array] = prop_SO(psi0, field, eH, mu_SO, tau, opt)
% Propagates states using the Split-Operator method.
arguments
	psi0  double {mustBeFinite}
	field double {mustBeFinite}
	eH    double {validateattributes(eH,{'double'},{'finite','square'})}
	mu_SO    {mustBeA(mu_SO,["struct","cell","numeric"])}
	tau   (1,1) double {mustBeFinite}
	opt.direction string {mustBeMember(opt.direction,["fwd","bwd"])} = "fwd";
end
[Npts, ~] = size(eH);
[Ndim, Ntau] = size(field);
[Nrow, Ncol] = size(psi0);

if (opt.direction=="fwd") & Nrow~=Npts
	% check if the initial state has valid dimensions
	error('ERROR: the number of rows in the initial state do not match the Hamiltonian')
elseif (opt.direction=="bwd") & Ncol~=Npts
	% test if the initial state has valid dimensions
	error('ERROR: the number of columns in the initial state do not match the Hamiltonian')
end

% check the diagonalization of mu, and compute if necessary
mu_SO = prop_SO_init(mu_SO);

if (size(mu_SO.U,3)~=Ndim+1) | (size(mu_SO.D,2)~=Ndim)
	% check if Mu has the correct number of components
	error("ERROR: the number of field components does not match the number of dipole moments")
elseif (size(mu_SO.D,1)~=Npts) | all(size(mu_SO.U,1)~=[1,Npts]) | all(size(mu_SO.U,2)~=[1,Npts])
	% check if Mu has the correct size
	error("ERROR: the size of the diagonalization does not match the system size")
end

switch opt.direction
	case "fwd"
		% initial states are column vectors
		psi_array = complex(zeros(Nrow,Ncol,Ntau+1));
		psi_array(:,:,1) = psi0;
		psi_array = prop_SO_fwd(psi_array, field, eH, mu_SO, tau);
	case "bwd"
		% initial states are row vectors
		psi_array = complex(zeros(Nrow,Ncol,Ntau+1));
		psi_array(:,:,Ntau+1) = psi0;
		psi_array = prop_SO_bwd(psi_array, field, eH, mu_SO, tau);
end
% done
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [psi_array] = prop_SO_fwd(psi_array, field, eH, mu_SO, tau)
% forwards-propagate a set of states using the Split-Operator method
[Ndim, Ntau] = size(field);
psi = psi_array(:,:,1);
orig_norm = sum((abs(psi).^2),1);

% reconfigure Mu.U for faster memory access
U = mat2cell(complex(mu_SO.U),size(mu_SO.U,1),size(mu_SO.U,2),ones(1,size(mu_SO.U,3)));
U = cellfun(@complex,U,"UniformOutput",false);

U{1} = eH * U{1};
U{Ndim+1} = U{Ndim+1} * eH;

for jj = 1:Ntau
	% time-dependent potential
	eD = exp(-1i * tau * field(:,jj).' .* mu_SO.D);
	% propagate	
	for kk = Ndim:-1:1
		psi = eD(:,kk) .* (U{kk+1} * psi);
		if kk==1; psi = U{kk} * psi; end
	end
	% normalize and save
	curr_norm = sum((abs(psi).^2),1);
	psi = psi .* sqrt(orig_norm ./ curr_norm);
	psi_array(:,:,jj+1) = psi;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [psi_array] = prop_SO_bwd(psi_array, field, eH, mu_SO, tau)
% backwards-propagate a set of states using the Split-Operator method
[Ndim, Ntau] = size(field);
psi = psi_array(:,:,Ntau+1);
orig_norm = sum((abs(psi).^2),2);

% reconfigure Mu.U for faster memory access
U = mat2cell(complex(mu_SO.U),size(mu_SO.U,1),size(mu_SO.U,2),ones(1,size(mu_SO.U,3)));
U = cellfun(@complex,U,"UniformOutput",false);

U{1} = eH * U{1};
U{Ndim+1} = U{Ndim+1} * eH;
U = cellfun(@ctranspose,U); 

for jj = flip(1:Ntau)
	% time-dependent potential
	eD = exp(-1i * tau * field(:,jj).' .* mu_SO.D);
	% propagate
	for kk = Ndim:-1:1
		psi = (psi * U{kk+1}) .* eD(:,kk)';
		if kk==1; psi = psi * U{kk}; end
	end
	% normalize and save
	curr_norm = sum((abs(psi).^2),2);
	psi = psi .* sqrt(orig_norm ./ curr_norm);
	psi_array(:,:,jj) = psi;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%