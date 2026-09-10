   __    ___          _____  ___    ___    ____        
 /‾_‾\ |‾   ‾||‾| |‾|(‾   ‾)|‾ _‾\ |‾ _‾\ |‾ __)|‾|    
( (_\_)| |‾| || | | | ‾| |‾ ( (_) )( (_) )| (_  | |    
 \__ \ | | | || | | |  | |  |    / |    / |  _) | |    
(‾\_) )|  (‾\)( |_| ) _| |_ | |\ \ | |\ \ | (__ | |__  
 \___/  \__\_) \___/ (_____)|_| \_)|_| \_)|____)|____) 

(Streamlined Quantum Unified Interface for Researching Real-time Excitations with Light)
------------------------------------------------------------------------------------------------------
The SQUIRREL software package was built in MATLAB R2023b to enable modeling of time-dependent quantum
effects in nano/mesoscale systems. As the geometry of these systems can be fairly complex, the 
SQUIRREL software package utilizes a finite-element approach to discretize the wavefunctions and other
spatially-dependent properties of the system. 

access the corresponding paper: https://doi.org/10.1016/j.cpc.2025.109861

BEFORE RUNNING
----------------------------------
Make sure the "functions" folder and its contents have been added to the MATLAB path (right-click
the "functions" folder and select: Add to Path > Selected Folders and Subfolders). 

Additionally, installation of the Partial Differential Equation Toolbox is required before calling 
the functions "meshify_geom.m" and "meshify_STL.m".

For more details on how to use the SQUIRREL package, one can access the included 
documentation located at docs/index.html, or refer to the included example scripts:

TUTORIAL FILES
----------------------------------
- example_meshing:		Demonstrates the process of generating meshes for 2D and 3D systems
- example_2D_nanowire:		Demonstrates time-independent and time-dependent calculations on the 
				2D cross-section of a core-shell nanowire
- example_anistropic_nanodot:	Sets up calculations for a 3D nanodot with anisotropic effective mass
- example_field_operators:	Visualization of the action of the electric/magnetic dipole moment on the 
				eigenstates of a 3D geometry
- example_multi_electron:	Demonstrates the setup of time-dependent calculations with interacting electrons

 Workflow:
----------------------------------
- define the geometry of the system / generate a mesh
- define spatially-dependent potentials/properties
- construct the Hamiltonian of the system
- use element-dropping to sparsify the Hamiltonian (if needed)
- compute eigenstates / other time-independent quantities
- define any time-dependent fields
- propagate the chosen wavefunction(s)
- display/visualize the results


LIST OF KEY FUNCTIONS: 
------------------------------------------------------------------------------------------------------

To generate a mesh:
----------------------------------
meshify_geom.m	Generates a 2D mesh from a geometry description with the requested number of nodes
meshify_STL.m	Generates a 2D/3D mesh from an STL file

To restructure/improve the quality of a mesh:
----------------------------------
mesh_anneal.m	Conditions 2D meshes by swapping edges within the triangulation
mesh_smooth.m	Conditions a mesh by averaging the location of nodes and their neighbors
mesh_symrcm.m	Reorders the mesh using symrcm to reduce the bandwidth of the adjacency matrix
mesh_symamd.m	Reorders the mesh using symamd to minimize the degree of the adjacency matrix


Constructing the FE system:
----------------------------------
atomic_units.m	Converts common SI units to their equivalent in Hartree atomic units
coeff_eval.m	Takes a function/string/double input and converts to a valid coeff array
FEmat.m		Constructs nullspace, overlap, orthobasis, vector/scalar potential, and stiffness matrices
FE_dipole.m	Returns the x,y,z components of the dipole-moment matrices for constant electric/magnetic fields
FE_init.m	Returns the main components of the Hamiltonian in the sparse/orthonormal basis
to_basis.m	Converts operators and states between the full <--> orthonormal or the full <--> sparse basis
eigenstates.m	Computes the requested number of lowest-energy eigenstates/values


Element-dropping / Sparsification methods
----------------------------------
drop_neighbor.m		Nearest-neighbor-based element dropping
drop_distance.m		Distance-based element dropping
drop_perturb.m		Perturbation-based element dropping
sparsify.m		Converts matrix to sparse if it improves memory and compute speed


Propagation Methods:
----------------------------------
prop_PA.m	Pade approximant propagator
prop_MH.m	Al-Mohy Higham propagator
prop_CN.m	Crank-Nicolson propagator (include the overlap matrix "M" for sparse propagation)
prop_SO_init.m	Initialize the Split-Operator propagator
prop_SO.m	Split-Operator propagator


Multi-Electron Systems:
----------------------------------
v_xc_LDA.m		Computes the exchange-correlation potential coefficients using the local-density approximation
setup_poisson.m		Sets up solver for the Poisson equation, needed to incorporate Coulombic interactions
setup_operator.m	Sets up a function that converts a potential into its operator form
prop_CN_multi.m		Crank-Nicolson propagator for multi-electron systems

Plotting:
----------------------------------
meshplot2D.m	Plots several sets of mesh data in the same figure; use the mouse scroll wheel 
		or gcf().UserData.setIndex() to change the displayed data
meshplot3D.m	meshplot2D for 3D systems, has several different display styles available

