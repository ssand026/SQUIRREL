function [out] = v_xc_LDA(mass,perm,rho,opt)
% Compute the exchange-correlation potential using the local-density approximation. 
% If the inputs are left empty, return a function for generating the exchange-correlation
% potential: out = func(mass,perm,dens). If the inputs "mass", "perm", and "dens"
%  are given numeric values, return the numeric exchange-correlation potential.
arguments
	mass (:,1) double {mustBeFinite,mustBeReal,mustBeNonzero} = [];
	perm (:,1) double {mustBeFinite,mustBeReal,mustBeNonzero} = [];
	rho (:,1) double {mustBeFinite,mustBeReal} = [];
	opt.method {mustBeMember(opt.method,[1,2])} = 2;
end

% There are two sets of empirical values for this formulation of the LDA exchange-correlation 
% energy, found in units of Rydberg in the papers below
switch opt.method
	case 1
		% L. Hedin and B. I. Lundqvist, "Explicit local exchange-correlation potentials,"
		% J. Phys. C: Solid State Phys., vol. 4, no. 14, p. 2064, Oct. 1971, 
		% doi: 10.1088/0022-3719/4/14/022.
		A = 0.0368;
		B = 21.0;
	case 2
		% O. Gunnarsson and B. I. Lundqvist, "Exchange and correlation in atoms, molecules, and 
		% solids by the spin-density-functional formalism,"
		% Phys. Rev. B, vol. 13, no. 10, pp. 4274–4298, May 1976,
		% doi: 10.1103/PhysRevB.13.4274
		A = 0.0545;
		B = 11.4;
end

% check that inputs are non-empty and of equal size (or scalars)
len = [numel(perm),numel(mass),numel(rho)];
doNumeric = all(len~=0) && all(len==1 | len==max(len));

if doNumeric
	% The quantity rs is the radius of the sphere (in bohr) that contains exactly one electron.
	% As the density can be zero-valued, we should avoid taking its inverse. Therefore we replace
	% all instances of 1/rs with inv_rs, which is proportional to the cube root of the density.
	inv_rs = (4*pi/3 * rho).^(1/3);
	
	% unitless scaling factor
	alpha = (4/(9*pi))^(1/3);

	% return the explicit potential (in Hartree)
	out = -2/(pi*alpha) * (inv_rs + A*log(1 + B*inv_rs)) .* (mass./(perm.^2));
else
	% unitless scaling factor
	alpha = (4/(9*pi))^(1/3);

	% return the potential generating function (in Hartree)
	out = @(mass,perm,rho) -2/(pi*alpha) .* (mass./(perm.^2)) .* ...
		( (4*pi/3*rho).^(1/3) + A*log( 1 + B*((4*pi/3*rho).^(1/3)) ) );
end
end
