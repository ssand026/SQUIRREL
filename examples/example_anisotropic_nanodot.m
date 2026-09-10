% Perform calculations of electron/electron-hole states with an anisotropic 
% effective mass in a Wurtzite GaN nanodot

% units
%---------------------------------------
fs = 1e-15 * atomic_units("seconds");
nm = 1e-9 * atomic_units("meters");
eV = 1 * atomic_units("eV");
mT = 1e-3 * atomic_units("tesla");

% create mesh with
[p,t] = meshify_STL("WZ_nanodot.stl",20000);
[p,t] = mesh_symrcm(p,t);
[p,t] = mesh_smooth(p,t,20);

p = 12.5*nm*p; % scale mesh

% input parameters
%---------------------------------------
q = -1; % charge

% the band-structure and other potentials must be scalar quantities
bandCoeff = @(z) (4 - (z-min(z))./(max(z)-min(z)))*eV;
bandCoeff = coeff_eval(p,t,bandCoeff);
meshplot3D(p,t, bandCoeff)

% however, material properties like the inverse of the effective mass or the 
% relative permittivity can be tensors to allow for anisotropic effects

massCoeff = [0.3, 0.2, 0.2]; % eff. mass along x,y,z

% or in tensor form, 1/2*massCoeff
invMassC = [1.67, 0.00, 0.00;
			0.00, 2.50, 0.00;
			0.00, 0.00, 2.50];
% the anisotropy can have spatial dependence too!

% construct matrices
%---------------------------------------
% Dirichlet nullspace matrix
B = FEmat(p,t,"dirichlet");

% overlap matrix
M = FEmat(p,t,"overlap");
M = to_basis("sparse",M,"mat",B);

% static potential matrix
V = FEmat(p,t,"potential",bandCoeff);
V = to_basis("sparse",V,"mat",B);

% kinetic energy matrix
K = FEmat(p,t,"stiffness",invMassC);
K = to_basis("sparse",K,"mat",B);

% static Hamiltonian
%---------------------------------------
H0 = K + V;
[states, energy, ~] = eigenstates(H0,M,"numEig",16);

% compute dipole moment matrices
%---------------------------------------
% compute electric field moments
muE = FE_dipole(p,t,"electric","charge",q);
for kk = 1:length(muE)
	muE{kk} = to_basis("sparse",muE{kk},"mat",B);
end

% compute magnetic field moments
[muB,muZ] = FE_dipole(p,t,"magnetic",invMassC,"charge",q);
for kk = 1:3
	muB{kk} = to_basis("sparse",muB{kk},"mat",B);
	muZ{kk} = to_basis("sparse",muZ{kk},"mat",B);
end

% plot eigenstates
%---------------------------------------
fig1 = meshplot3D(p,t,B*states,"plotStyle","surface");
view(0,90);
drawnow


%% propagate under an electromagnetic field
%---------------------------------------
Ntau = 800;
tau = 0.1*fs;
psi0 = states(:,11);

% setup the time-dependent electromagnetic field
Ex = 5e-2 * eV/nm .* sinpi(linspace(0,4,Ntau));
Ey = 3e-2 * eV/nm .* cospi(linspace(0,4,Ntau));
Bz = 5e2 * mT .* sinpi(linspace(0,1,Ntau)).^2;

% propagate (WARNING: MAY TAKE A FEW MINUTES TO COMPLETE)
field = [Ex; Ey; Bz; Bz.^2];
mu = {muE{1}, muE{2}, muB{3}, muZ{3}};
[psi_prop] = prop_CN(psi0, field, H0, mu, tau, M);

% plot time-evolution
rhoT = abs(squeeze(psi_prop)).^2;
fig2 = meshplot3D(p,t,B*rhoT,"plotStyle","volume","fixedCLims",true);
view(-10,40);

% display animation of the time-evolution
%---------------------------------------
framet = 1/30;
for kk = 1:10:size(rhoT,2)
	tic;
	fig2.UserData.setIndex(kk);
	drawnow;
	tim = toc;
	if tim<framet
		pause(framet-tim)
	end
end

