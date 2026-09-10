function [K,V,M,B,muE,muB,muZ] = FE_init(pts,tri,massCoeff,bandCoeff,basis)
% Computes all components of the finite-element Hamiltonian in the specified basis.  
% The input coefficients "massCoeff" and "bandCoeff" can be defined using a
% function_handle/string/char, or via a numeric array with valid dimensions.
arguments
	pts (:,:) {mustBeFinite,mustBeReal}
	tri (:,:) {mustBeInteger,mustBePositive}
	massCoeff {mustBeA(massCoeff,{'function_handle','char','string','numeric','logical'})}
	bandCoeff {mustBeA(bandCoeff,{'function_handle','char','string','numeric','logical'})}
	basis string {mustBeMember(basis,["sparse","orth"])}
end

switch basis
	case "orth"
		% orthonormal basis conversion matrix
		B = FEmat(pts,tri,"orthobasis");

		% kinetic energy matrix
		K = FEmat(pts,tri,"stiffness",massCoeff);
		K = to_basis("orth",full(K),"mat",B);

		% static potential matrix
		V = FEmat(pts,tri,"potential",bandCoeff);
		V = to_basis("orth",full(V),"mat",B);

		% add the overlap/full-space-conversion matrices
		M = eye(size(V),'like',V);

		% compute electric field moments
		if nargout>=5
			muE = FE_dipole(pts,tri,"electric");
			muE = cellfun(@(x) to_basis("orth",full(x),"mat",B), muE, "UniformOutput",false);
		end

		% compute the magnetic field moments
		if nargout>=6
			[muB,muZ] = FE_dipole(pts,tri,"magnetic",massCoeff);
			muB = cellfun(@(x) to_basis("orth",full(x),"mat",B), muB, "UniformOutput",false);
			muZ = cellfun(@(x) to_basis("orth",full(x),"mat",B), muZ, "UniformOutput",false);
		end

	case "sparse"
		% dirichlet nullspace matrix
		B = FEmat(pts,tri,"dirichlet");

		% kinetic energy matrix
		K = FEmat(pts,tri,"stiffness",massCoeff);
		K = to_basis("sparse",K,"mat",B);

		% static potential matrix
		V = FEmat(pts,tri,"potential",bandCoeff);
		V = to_basis("sparse",V,"mat",B);

		% overlap matrix
		M = FEmat(pts,tri,"overlap");
		M = to_basis("sparse",M,"mat",B);

		% compute the electric field moments
		if nargout>=5
			muE = FE_dipole(pts,tri,"electric");
			muE = cellfun(@(x) to_basis("sparse",x,"mat",B), muE, "UniformOutput",false);
		end

		% compute the magnetic field moments
		if nargout>=6
			[muB,muZ] = FE_dipole(pts,tri,"magnetic",massCoeff);
			muB = cellfun(@(x) to_basis("sparse",x,"mat",B), muB, "UniformOutput",false);
			muZ = cellfun(@(x) to_basis("sparse",x,"mat",B), muZ, "UniformOutput",false);
		end
end