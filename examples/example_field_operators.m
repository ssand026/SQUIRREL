% EXAMPLE_FIELD_OPERATORS
% displays the action of the z-oriented field operators on the eigenstates of a 3D system

% units
%---------------------------------------
fs = 1e-15 * atomic_units("seconds");
nm = 1e-9 * atomic_units("meters");
eV = 1 * atomic_units("eV");
mT = 1e-3 * atomic_units("tesla");

% load in a cone-shaped geometry
[p,t] = meshify_STL("cone.stl",10000);
[p,t] = mesh_smooth(p,t,30);
p = 20*nm*p;

% compute system matrices
%---------------------------------------
massCoeff = 0.2;
bandCoeff = 0;

[K,V,M,B,muE,muB,muZ] = FE_init(p,t,1./(2*massCoeff),bandCoeff,"sparse");
H0 = K + V;

% compute eigenstates and plot
%---------------------------------------
[states, energy, didConv] = eigenstates(H0,M,"numEig",32);
fig1 = meshplot3D(p,t,B*states,"plotStyle","surface"); view(-10,40);

% compute the action of the eigenstates on the dipole moments
%---------------------------------------
dim = 3;
psiE = M\(muE{dim} * states);
psiB = M\(muB{dim} * states); psiB = imag(psiB);
psiZ = M\(muZ{dim} * states);

% plot the action of the field operators on the eigenstates
fig2 = meshplot3D(p,t,B*psiE,"plotStyle","surface","cmapLimits",["-2*std","+2*std"],"fixedCLims",true);
fig3 = meshplot3D(p,t,B*psiB,"plotStyle","surface","cmapLimits",["-2*std","+2*std"],"fixedCLims",true);
fig4 = meshplot3D(p,t,B*psiZ,"plotStyle","surface","cmapLimits",["-2*std","+2*std"],"fixedCLims",true);
drawnow


%% propagate under a magnetic field
%---------------------------------------
Ntau = 400;
tau = 0.5*fs;
psi0 = states(:,3);

% define the time-dependent magnetic field
Bz = 25e3*mT .* (sinpi(linspace(0,5,Ntau)).^2);

% propagate (WARNING: MAY TAKE A FEW MINUTES TO COMPLETE)
[psi_prop] = prop_CN(psi0,[Bz; Bz.^2],H0,{muB{3},muZ{3}},tau,M);

% plot time-evolution
%---------------------------------------
rhoT = to_basis("full",squeeze(abs(psi_prop).^2),"col",B);
meshplot3D(p,t,rhoT,"plotStyle","surface","cmapLimits",["0","max"]);

