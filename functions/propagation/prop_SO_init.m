function [mu_SO] = prop_SO_init(mu, diagMethod)
% Initialization process for the Split-Operator propagator.
%
% Diagonalizes the time-dependent components of the Hamiltonian for fast 
% exponentiation in the split-operator process. For the dipole matrix mu, returns 
% mu_SO such that expm(b*mu) = mu_SO.U(:,:,1) * exp(b*mu_SO.D) * mu_SO.U(:,:,2)
arguments
mu {mustBeA(mu,["struct","cell","numeric"])}
diagMethod string {mustBeMember(diagMethod,["full","lazy","rowsum","colsum"])} = "full";
end

% parse input types
if isstruct(mu)
	% check if struct input has a valid form
	if ~(isfield(mu,"U") & isfield(mu,"D"))
		error("ERROR: struct does not contain the requisite fields")
	elseif size(mu.D,2)+1 ~= size(mu.U,3)
		error("ERROR: fields of Mu have diferent dimensonality")
	elseif any(size(mu.D,1)~=size(mu.U,[1 2])) & all(size(mu.U,[1 2])~=1)
		error("ERROR: fields of Mu have diferent sizes")
	else
		mu_SO = mu;
		return;
	end
elseif iscell(mu)
	Ndim = length(mu);    % number of dipole matrices
	Npts = length(mu{1}); % number of matrix elements
	% check size of input
	for kk = 1:Ndim
		if any(size(mu{kk})~=Npts)
			error("ERROR: the cells of mu must be square matrices of equal size")
		end
	end
elseif isnumeric(mu)
	Ndim = size(mu,3); % number of dipole matrices
	Npts = size(mu,1); % number of matrix elements
	% check size of input
	if size(mu,1)~=size(mu,2) | Ndim > 3
		error("ERROR: mu must be a square matrix or pagearray")
	end
	mu = mat2cell(mu,1,1,ones(1,Ndim));
else
	error("ERROR: mu must be a struct, cell, pagearray, or matrix")
end

% construct diagonalization
muDiags = zeros(Npts,Ndim);
muPerms = repmat([1],1,1,Ndim+1);

switch diagMethod
	case "full"
		% compute full diagonalization
		muPerms = repmat(eye(Npts),1,1,Ndim+1);
		for kk = 1:Ndim
			[V,D] = eig(double(mu{kk}),"vector");
			muDiags(:,kk) = reshape(D,[],1);
			
			muPerms(:,:,kk) = muPerms(:,:,kk) * V;
			muPerms(:,:,kk+1) = V \ muPerms(:,:,kk+1);
		end
	case "lazy"
		% extract diagonal component
		for kk = 1:Ndim
			D = diag(mu{kk});
			muDiags(:,kk) = reshape(D,[],1);
		end
	case "rowsum"
		% use row sum
		for kk = 1:Ndim
			D = sum(mu{kk},1);
			muDiags(:,kk) = reshape(D,[],1);
		end
	case "colsum"
		% use column sum
		for kk = 1:Ndim
			D = sum(mu{kk},1);
			muDiags(:,kk) = reshape(D,[],1);
		end
end

% add the diagonal/eigenvector matrices to the output
mu_SO.U = muPerms;
mu_SO.D = muDiags;
end