function [faces] = mesh_faces(pts,tri)
% Returns the list of all unique faces contained within the mesh
arguments
	pts (:,:) double {mustBeFinite,mustBeReal}
	tri (:,:) {mustBeInteger,mustBePositive}
end
Ndim = size(pts,1);
Ntri = size(tri,2);
facePerms = ncombsk(Ndim+1,Ndim);
Nperm = size(facePerms,1);

faces = zeros(Ndim, Nperm*Ntri);
for nn = 1:Nperm
	indx = (nn-1)*Ntri + [1:Ntri];
	faces(:,indx) = tri(facePerms(nn,:), :);
end
faces = unique(sort(faces,1).',"rows");
end