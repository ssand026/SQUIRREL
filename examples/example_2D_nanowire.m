% COMPUTE EIGENSTATES AND COMPARE PROPAGATORS FOR A 2D NANOWIRE CROSS-SECTION

% units
%---------------------------------------
fs = 1e-15 * atomic_units("seconds");
nm = 1e-9 * atomic_units("meters");
meV = 1e-3 * atomic_units("eV");

% choose a ground-state configuration to load
% (results obtained using the HADOKEN software package)
%---------------------------------------
load("hex_5k.mat")
% load("triG_3k.mat")
% load("triN_3k.mat")

% compute system matrices
[K,V,M,B,muE] = FE_init(p,t,1./(2*massCoeff),Vcb+Vee+Vxc,"orth");
H0 = K + V;

% compute eigenstates and plot
%---------------------------------------
[states,energy] = eigenstates(H0,"numEigs",16);
states_full = to_basis("full",states,"col",B);

meshplot2D(p,t,states_full,"cmapLimits",["-mag","+mag"]);

%% propagate
%---------------------------------------
% define E-field
Ntau = 500;
tau = 1*fs;
Exy = [sinpi(2*linspace(0,5,Ntau)); cospi(2*linspace(0,5,Ntau))];
Exy = 1*meV/nm * Exy;

psi0 = states(:,4); % choose inital state

% SO-propagation
mu_SO = prop_SO_init(muE);
eH = expm(-1i*tau/2*H0);
[psi_SO] = prop_SO(psi0, Exy, eH, mu_SO, tau);

% drop elements from the Hamiltonian for CN method
H0_drop = sparse(drop_perturb(H0,1e-12));
for ii = 1:numel(muE)
	muE_drop{ii} = sparse(drop_perturb(muE{ii},1e-12)); 
end
fprintf("Dropped Hamiltonian density: %.3f"+"\n",nnz(H0_drop)/numel(H0));

% CN-propagation (WARNING: WILL TAKE SEVERAL MINUTES TO COMPLETE)
[psi_CN] = prop_CN(psi0, Exy, H0_drop, muE_drop, tau);

% compute probability density
rho_SO = abs(to_basis("full",psi_SO,"col",B)).^2;
rho_CN = abs(to_basis("full",psi_CN,"col",B)).^2;

% display time-dependent probability density
meshplot2D(p,t,rho_SO,"cmapLimits",["eps","avg+3*std"],"fixedCLims",true);

% display the difference between SO/CN propagation
meshplot2D(p,t,rho_SO-rho_CN,"cmapLimits",["-4*std","+4*std"],"fixedCLims",false);