function [a1,a2,a3] = mesh_angles(pts,tri)
% Returns the angles for each triangle in the mesh
arguments
	pts (2,:) double {mustBeFinite,mustBeReal}
	tri {mustBeInteger,mustBePositive}
end
vec2angle = @(a,b) acos(sum(a.*b,2)./(sqrt(sum(a.^2,2)).*sqrt(sum(b.^2,2))));

n1 = pts(:,tri(1,:)).';
n2 = pts(:,tri(2,:)).';
n3 = pts(:,tri(3,:)).';

a1 = vec2angle(n2-n1,n3-n1);
a2 = vec2angle(n1-n2,n3-n2);
a3 = vec2angle(n1-n3,n2-n3);
end