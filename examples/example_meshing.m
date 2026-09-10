% EXAMPLE_MESHING

% The first step in any calculation is to generate a mesh for the system of 
% interest. Often it is preferable to define a geometry for your system and have 
% the mesh generated automatically.

% Here we use the parameter num_nodes to control the size of the mesh. While a 
% larger number of mesh nodes can improve the accuracy of the simulation, the 
% size of the basis (and therefore the size of the matrices) scales proportionally,
% resulting in higher computational costs. For non-sparse calculations, we do not
% recommend setting num_nodes above ~5000 for general use. For fully-sparse 
% calculations, one can easily use meshes with over 10,000 nodes.

num_nodes = 2000;

% The Partial Differential Equation Toolbox has a couple of built-in tools for defining 
% custom geometries, such as the function "decsg.m" or the PDE Modler app. In 
% addition, the PDE toolbox includes a variety of prebuilt geometries, which can 
% be called via: @lshapeg, @cardg, @circleg, @cirsg, @crackg, @scatterg, or @squareg.
%
% Currently, these approaches only support the generation of 2D geometries. 

% generate a mesh using a built-in geometry
geom = @lshapeg;
[p,t] = meshify_geom(geom,num_nodes);
meshplot2D(p,t,[])

% define a custom geometry using decsg.m
c1 = [1;          0;    1; 0.9]; % circle 1
c2 = [1; +sqrt(3)/2; -1/2; 0.9]; % circle 2
c3 = [1; -sqrt(3)/2; -1/2; 0.9]; % circle 3

geom = decsg([c1,c2,c3],'(A+B+C)',['ABC']);
geom = csgdel(geom,1);
[p,t] = meshify_geom(geom,num_nodes);
meshplot2D(p,t,[])

% to see the quality of a given mesh, one can use the command:
meshplot2D(p,t,pdetriq(p,t),"cmapLimits",["min","1"])

% For 3D geometries, the SQUIRREL function meshify_STL.m, allows one to use a
% STL file to create a geometry. The resulting mesh will have a maximum length
% of one across all axes.
[p,t] = meshify_STL("cone.stl",num_nodes);
meshplot3D(p,t,[])