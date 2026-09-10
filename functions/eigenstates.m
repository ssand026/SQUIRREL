function [states, energy, didConv] = eigenstates(H,M,opt)
% Returns the normalized eigenvectors and the eigenvalues for a Hamiltonian system.
% Specifying the overlap-matrix "M" as the 2nd argument solves for the eigen-
% values/vectors H*v = λ*M*v instead of H*v = λ*v
arguments
	H (:,:) double {validateattributes(H,{'double'},{'finite','square'})}
	M (:,:) double {validateattributes(M,{'double'},{'finite','square'})} = [];
	opt.numEigs (1,1) double {mustBeInteger,mustBePositive} = 32;
end
num_eigs = min(size(H,1),opt.numEigs);

if isempty(M)
	% orthonormal basis
	[states,energy,flag] = eigs(H,num_eigs,"smallestreal","FailureTreatment","keep",...
		"Tolerance",1e-12,"SubspaceDimension",min(max(3*num_eigs,32),size(H,1)));
	states = states ./ sqrt(sum(abs(states).^2,1));
else
	% sparse basis
	[states,energy,flag] = eigs(H,M,num_eigs,"smallestreal","FailureTreatment","keep", ...
		"Tolerance",1e-12,"SubspaceDimension",min(max(3*num_eigs,32),size(H,1)));
	states = states ./ sqrt(sum(conj(states) .* (M * states),1));
end
energy = reshape(diag(energy),[],1);
didConv = ~flag;
if any(flag); warning('some eigenstates did not converge'); end
end