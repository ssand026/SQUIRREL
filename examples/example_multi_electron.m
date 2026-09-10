% MULTI-ELECTRON TUTORIAL
% Showcase propagation of multi-electron systems

% define units
nm = 1e-9 * atomic_units("meters");
eV = atomic_units("eV");
fs = 1e-15 * atomic_units("seconds");

% choose a ground-state configuration to load
% (results obtained using the HADOKEN software package)
load("hex_5k.mat")
% load("triG_3k.mat")
% load("triN_3k.mat")

% plot the occupied eigenstates
meshplot2D(p,t,psi,"cmapLimits",["-mag","+mag"],"fixedCLims",true); colorbar;
drawnow;

% compute components of the Hamiltonian
basis = "sparse";
[K,V_CB,M,B,muE,muB,muZ] = FE_init(p,t,1./(2*massCoeff),bandCoeff,"sparse");

% initialize Poisson solver
q = -1; % electron charge
poissFunc = setup_poisson(p,t,permCoeff,"method","decompose");
to_operator = setup_operator(p,t);
to_full = @(data,shape) to_basis("full",data,shape,B);
to_redu = @(data,shape) to_basis("sparse",data,shape,B);

% compute the static Coulombic potential
n0 = coeff_check(p,t,rho_0); % 
nd = coeff_check(p,t,rho_pos); % donor d
v0 = q*poissFunc("rho",n0+nd); % coefficient form
V0 = to_redu(to_operator(v0),"mat"); % operator form

% setup the dynamic coulombic potentials
n_e = @(ket) q * sqrt(2*massCoeff)/pi .* (abs(ket).^2 * sqrt(E_fermi-E_schrod));
vee = @(ket) q*poissFunc("rho",n_e(ket));

% combine into the time-dependent Kohn-Sham potential
vxc = @(ket) q^2*v_xc_LDA(massCoeff,permCoeff,sum(abs(ket).^2,2));
vks = @(ket) vee(to_full(ket,"col")) + vxc(to_full(ket,"col")); % coefficient form
Vks = @(ket) to_redu(to_operator(vks(ket)),"mat"); % operator form

% Define the time-independent component of the Hamiltonian
H0 = K + V_CB + V0;

% plot the total static potential
meshplot2D(p,t,bandCoeff+v0,"cmapLimits",["min","max"]); colorbar
drawnow;

% plot the initial Kohn-Sham potential
meshplot2D(p,t,vks(psi),"cmapLimits",["min","max"]); colorbar
drawnow;

%% propagate

% setup field and timesteps
tau = 0.1*fs;
Ntau = 600;
field = [0.03; 0.02] .* eV/nm .* sinpi([4;5].*linspace(0,1,Ntau));

% run propagation
psi0 = to_basis("sparse",psi,"col",B);
psi_t = prop_CN_multi(psi0, field, H0, muE, tau, M, Vks);

% convert psi_t to the full-basis probability density
psi_t = pagemtimes(full(B),psi_t);
rho_t = squeeze(sum(abs(psi_t).^2,2));

% plot
meshplot2D(p,t,rho_t,"cmapLimits",["0","mag"]); colorbar;
drawnow;


