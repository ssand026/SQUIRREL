function [varargout] = FE_dipole(pts,tri,fieldType,massCoeff,opt)
% Constructs the dipole-moment matrices for spatially-invariant electric or magnetic fields.
% The output will be a cell array, where each entry corresponds to the dipole
% operator for the field along each axis.
arguments
	pts (:,:) {mustBeFinite,mustBeReal}
	tri (:,:) {mustBeInteger,mustBePositive}
	fieldType string {mustBeMember(fieldType,["electric","magnetic"])}
	massCoeff = 1/2;
	opt.charge (1,1) double {mustBeNumeric,mustBeReal,mustBeFinite} = -1;
end
[Ndim,Npts] = size(pts);
q = opt.charge;

switch fieldType
	case "electric"
		% electric dipole moments
		muE = cell(1,Ndim);
		for dd = 1:Ndim
			muE{dd} = q * FEmat(pts,tri,"potential",-pts(dd,:));
		end
		varargout{1} = muE;

	case "magnetic"
		% magnetic dipole moments
		if Ndim~=2 && Ndim~=3
			error('ERROR: system must be 2D or 3D for magnetic moments')
		elseif Ndim==2
			% z-direction B
			vB = -1/2 * [pts(2,:); -pts(1,:)];
			muB{1} = 1i*q*FEmat(pts,tri,"skew-vector",massCoeff,vB);
			muZ{1} = q^2 *FEmat(pts,tri,"symm-vector",massCoeff,vB);
		elseif Ndim==3
			% x-direction B
			vB = -1/2 * [zeros(1,Npts); pts(3,:); -pts(2,:)];
			muB{1} = 1i*q*FEmat(pts,tri,"skew-vector",massCoeff,vB);
			muZ{1} = q^2 *FEmat(pts,tri,"symm-vector",massCoeff,vB);

			% y-direction B
			vB = -1/2 * [-pts(3,:); zeros(1,Npts); pts(1,:)];
			muB{2} = 1i*q*FEmat(pts,tri,"skew-vector",massCoeff,vB);
			muZ{2} = q^2 *FEmat(pts,tri,"symm-vector",massCoeff,vB);

			% z-direction B
			vB = -1/2 * [pts(2,:); -pts(1,:); zeros(1,Npts)];
			muB{3} = 1i*q*FEmat(pts,tri,"skew-vector",massCoeff,vB);
			muZ{3} = q^2 *FEmat(pts,tri,"symm-vector",massCoeff,vB);
		end
		varargout{1} = muB;
		varargout{2} = muZ;
end
% done
end