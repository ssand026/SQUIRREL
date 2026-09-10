function [bnodes] = boundary_nodes(pts,tri)
% Returns the logical index corresponding to nodes on the mesh boundary
arguments
	pts (:,:) double {mustBeFinite,mustBeReal}
	tri (:,:) {mustBeInteger,mustBePositive}
end
[Ndim, Npts] = size(pts);
dt = sort(tri(1:Ndim+1,:),1);

% set of face-node permutations
per = repmat([1:Ndim+1],1,Ndim+1);
per(1:Ndim+2:end) = [];
faces = reshape(dt(per,:),Ndim,[]);

% use minimal memory to store face indices
data_dim = Ndim*log2(Npts);
if data_dim < 8
	new_type = "uint8";
elseif data_dim < 16
	new_type = "uint16";
elseif data_dim < 32
	new_type = "uint32";
elseif data_dim < 64
	new_type = "uint64";
else
	% linearized indices would exceed memory limit (should be difficult to do)
	error('ERROR: number of mesh nodes is too large to linearize face indices')
end

% convert face indexes to a linear index
faces = cast(faces-1,new_type);
scl = cast([Npts.^[0:Ndim-1]],new_type);

faceID = scl(1) * faces(1,:);
for jj = 2:length(scl)
	faceID = faceID + (scl(jj) * faces(jj,:));
end

% find faces belonging to a single simplex
faceID = sort(faceID);
test = (faceID ~= [faceID(2:end),0]) & (faceID ~= [0,faceID(1:end-1)]);
faceID = faceID(test);

% convert faceID back to node indices
bnodes = false(Npts,1);
for kk = Ndim:-1:1
	r = mod(faceID,scl(kk));
	indx = (faceID - r)/scl(kk);
	indx = indx+1;
	bnodes(indx) = true;
	faceID = r;
end
end