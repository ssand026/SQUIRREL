function [out] = FEmat(pts, tri, matType, scalarC, vectorC)
% Returns the specified finite-element matrix
%
% INPUTS:
%	pts: mesh node locations
%	tri: mesh triangulation
%	matType (string): type of matrix to construct
%	scalarC: scalar/tensor coefficient for the matrix
%	vectorC: vector coefficient for the matrix
%
% SEE ALSO: FEMAT_NULLSPACE, FEMAT_ORTHOBASIS, FEMAT_OVERLAP, PARPOOL
arguments
	pts (:,:) double {mustBeFinite,mustBeReal}
	tri (:,:) {mustBeInteger,mustBePositive}
	matType string {mustBeMember(matType,["potential","stiffness","surface",...
		"skew-vector","symm-vector","overlap","dirichlet","nullspace","orthobasis"])}
	scalarC {mustBeA(scalarC,{'cell','string','char','function_handle','numeric','logical'})} = [];
	vectorC {mustBeA(vectorC,{'cell','string','char','function_handle','numeric','logical'})} = [];
end

% trim the triangulation to the nearest polynomial order
Ndim = size(pts,1);
vtx_len = cumsum(factorial(Ndim+1)./(factorial(Ndim-[0:Ndim]) .* factorial([1:Ndim+1])));
Npoly = find(size(tri,1)>=vtx_len,1,"last");
tri = tri(1:vtx_len(Npoly),:);

% evaluate the coefficients and check their size
switch matType
	case "potential"
		% eval scalar coeff
		[scalarC,order] = coeff_eval(pts,tri,scalarC);
		if order ~= 0
			error("ERROR: the scalar coefficient must be of order 0")
		end
	case "stiffness"
		% eval scalar/tensor coeff
		[scalarC,order] = coeff_eval(pts,tri,scalarC);
		if order >= 3
			error("ERROR: the scalar coefficient must be of order 0, 1, or 2")
		end
	case "surface"
		% eval scalar/tensor coeff
		[scalarC,order] = coeff_eval(pts,tri,scalarC);
		if order >= 3
			error("ERROR: the scalar coefficient must be of order 0, 1, or 2")
		end
		if ~isempty(vectorC)
			% check vector coeff
			[vectorC,order] = coeff_eval(pts,tri,vectorC);
			if order ~= 1
				error("ERROR: the vector coefficient must be of order 1")
			end
		end

	case {"skew-vector","symm-vector"}
		% eval scalar/tensor coeff
		[scalarC,order] = coeff_eval(pts,tri,scalarC);
		if order >= 3
			error("ERROR: the scalar coefficient must be of order 0, 1, or 2")
		end
		% check vector coeff
		[vectorC,order] = coeff_eval(pts,tri,vectorC);
		if order ~= 1
			error("ERROR: the vector coefficient must be of order 1")
		end
end

% compute appropriate matrix type
switch matType
	case "overlap"
		out = femat_overlap(pts,tri);
	case {"dirichlet","nullspace"}
		out = femat_nullspace(pts,tri);
	case "orthobasis"
		out = femat_orthobasis(pts,tri);
	case "potential"
		out = femat_potential(pts,tri,scalarC);
	case "stiffness"
		out = femat_stiffness(pts,tri,scalarC);
	case "surface"
		out = femat_surface(pts,tri,scalarC,vectorC);
	case "skew-vector"
		out = femat_skew_vector(pts,tri,scalarC,vectorC);
	case "symm-vector"
		out = femat_symm_vector(pts,tri,scalarC,vectorC);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [V] = femat_potential(pts,tri,scalarC)
% constructs the finite-element scalar-potential matrix
arguments
	pts (:,:) double {mustBeFinite,mustBeReal}
	tri (:,:) {mustBeInteger,mustBePositive}
	scalarC double {mustBeFinite}
end
[Ndim,Npts] = size(pts);
[Nvtx,Ntri] = size(tri);

% construct simplex coefficients
scalarC = coeff_expand(pts,tri,scalarC);

% compute integrand values
integrand = permute(scalarC,[3,4,1,2]);

% compute integration coefficients
svol = simplex_vol(pts,tri);
[intCoeff,ii,jj,kk] = bary_integral(Ndim,3);
ii_jj = sub2ind([Nvtx,Nvtx],ii,jj);

% loop over simplexes
[rows,cols] = npermsk(Nvtx,2);
rows = tri(rows(:),:);
cols = tri(cols(:),:);
vals = zeros(size(rows));

for ss = 1:Ntri
	int_ss = integrand(:,ss);
	Vij = accumarray(ii_jj,intCoeff .* int_ss(kk,:));
	vals(:,ss) = svol(ss) .* Vij(:);
end
V = sparse(rows,cols,vals,Npts,Npts);
V = (V + V.')/2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [K] = femat_stiffness(pts,tri,scalarC)
% constructs the finite-element stiffness matrix
arguments
	pts (:,:) double {mustBeFinite,mustBeReal}
	tri (:,:) {mustBeInteger,mustBePositive}
	scalarC double {mustBeFinite}
end
[Ndim,Npts] = size(pts);
[Nvtx,Ntri] = size(tri);

% construct simplex coefficients
baryvec = coeff_expand(pts,tri,baryvectors(pts,tri));
scalarC = coeff_expand(pts,tri,scalarC);

% compute integrand values
[aa,bb,cc] = npermsk(Nvtx,3);
v_aa = baryvec(:,:,aa,:);
m_bb = scalarC(:,:,bb,:);
v_cc = baryvec(:,:,cc,:);

if isequal(size(m_bb,[1,2]),[Ndim,Ndim])
	v_cc = permute(v_cc,[2,1,3,4]);
	integrand = sum(v_aa .* m_bb .* v_cc,[1,2]);
	integrand = permute(integrand,[3,4,1,2]);
else
	integrand = sum(v_aa .* m_bb .* v_cc,[1,2]);
	integrand = permute(integrand,[3,4,1,2]);
end

% compute integration coefficients
svol = simplex_vol(pts,tri);
[intCoeff,ii,jj,kk] = bary_integral(Ndim,3,[1 2]); % integration coeffs ignoring ii and jj
ii_jj = sub2ind([Nvtx,Nvtx],ii,jj);
ii_kk_jj = sub2ind([Nvtx,Nvtx,Nvtx],ii,kk,jj);

% loop over simplexes
[rows,cols] = npermsk(Nvtx,2);
rows = tri(rows(:),:);
cols = tri(cols(:),:);
vals = zeros(size(rows));

for ss = 1:Ntri
	int_ss = integrand(:,ss);
	Kij = accumarray(ii_jj, intCoeff .* int_ss(ii_kk_jj,:));
	vals(:,ss) = svol(ss) .* Kij(:);
end
K = sparse(rows,cols,vals,Npts,Npts);
K = (K + K.')/2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [B] = femat_surface(pts,tri,scalarC,vectorC)
% constructs the finite-element surface integral matrix
arguments
	pts (:,:) double {mustBeFinite,mustBeReal}
	tri (:,:) {mustBeInteger,mustBePositive}
	scalarC double {mustBeFinite}
	vectorC double {mustBeFinite} = [];
end
[Ndim,Npts] = size(pts);
svol = simplex_vol(pts,tri);

% construct simplex coefficients
bvec = coeff_expand(pts,tri,baryvectors(pts,tri));
scalarC = coeff_expand(pts,tri,scalarC);

% use baryvectors if vectorCoeff not given
if isempty(vectorC)
	vectorC = bvec;
	ignore = 2;
else
	vectorC = coeff_expand(pts,tri,vectorC);
	ignore = [];
end

% find simplices with faces on the domain boundary
isbndnode = boundary_nodes(pts,tri);
isbndtri = sum(isbndnode(tri(1:Ndim+1,:)),1)==Ndim;
bndtri = tri(1:Ndim+1,isbndtri);

% reduce the coefficients to only include bounding simplexes
svol = svol(:,isbndtri);
bvec = bvec(:,:,:,isbndtri);
scalarC = scalarC(:,:,:,isbndtri);
vectorC = vectorC(:,:,:,isbndtri);

% find boundary faces on a corner
iscorner = (accumarray(tri(:),ones(numel(tri),1),[Npts,1])==1);
iscornertri = any(iscorner(bndtri),1);
cfaces = mesh_faces(pts,bndtri(:,iscornertri));
cfaces = cfaces(any(iscorner(cfaces),2)&all(isbndnode(cfaces),2),:);

% find all other boundary faces
bfaces = mesh_faces(pts,bndtri(:,~iscornertri));
bfaces = bfaces(all(isbndnode(bfaces),2),:);

% join faces and extract relevant coefficient values
faces = [cfaces; bfaces].';
Nface = size(faces,2);

% find face areas/normals/coefficients
f_vol = zeros(1,Nface);
f_scl = zeros([size(scalarC,[1,2]),[Ndim,Nface]]);
f_vec = zeros([size(vectorC,[1,2]),[Ndim,Nface]]);
f_norm = zeros(Ndim,1,1,Nface);

for ff = 1:Nface
	% find which vertices correspond to the face
	isFaceVtx = ismember(bndtri,faces(:,ff));
	% find which simplex contains the specified face
	val = find(sum(isFaceVtx,1)==Ndim);
	if numel(val)~=1
		% the valid simplex must also contain a corner node
		val = val(any(iscorner(bndtri(:,val)),1));
	end
	% populate face parameters
	f_vol(:,ff) = svol(:,val);
	f_scl(:,:,:,ff) = scalarC(:,:,isFaceVtx(:,val),val);
	f_vec(:,:,:,ff) = vectorC(:,:,isFaceVtx(:,val),val);
	f_norm(:,:,:,ff) = -bvec(:,:,~isFaceVtx(:,val),val);
end
f_area = 2 .* squeeze(vecnorm(f_norm,2,1)) .* f_vol;
f_norm = f_norm ./ vecnorm(f_norm,2,1);

% compute integrand values
[aa,bb] = npermsk(Ndim,2);
if isequal(size(f_scl,[1,2]),[Ndim,Ndim])
	f_vec = permute(f_vec,[2,1,3,4]);
end
integrand = sum(f_norm .* f_scl(:,:,bb,:) .* f_vec(:,:,aa,:),[1,2]);
integrand = permute(integrand,[3,4,1,2]);

% compute integration coefficients
[intCoeff,ii,jj,kk] = bary_integral(Ndim-1,3,ignore);
ii_jj = sub2ind([Ndim,Ndim],ii,jj);
jj_kk = sub2ind([Ndim,Ndim],jj,kk);

% loop over faces
[rows,cols] = npermsk(Ndim,2);
rows = faces(rows,:);
cols = faces(cols,:);
vals = zeros(size(rows));

for ff = 1:Nface
	int_ff = integrand(:,ff);
	Bij = accumarray(ii_jj, intCoeff .* int_ff(jj_kk,:));
	vals(:,ff) = f_area(ff) .* Bij(:);
end
B = sparse(rows,cols,vals,Npts,Npts);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X] = femat_skew_vector(pts,tri,scalarC,vectorC)
% constructs the finite-element matrix for the complex-valued vector-potential operator (1st order)
arguments
	pts (:,:) double {mustBeFinite,mustBeReal}
	tri (:,:) {mustBeInteger,mustBePositive}
	scalarC double {mustBeFinite}
	vectorC double {mustBeFinite}
end
[Ndim,Npts] = size(pts);
[Nvtx,Ntri] = size(tri);

% construct simplex coefficients
baryvec = coeff_expand(pts,tri,baryvectors(pts,tri));
scalarC = coeff_expand(pts,tri,scalarC);
vectorC = coeff_expand(pts,tri,vectorC);

% compute integrand values
[aa,bb,cc] = npermsk(Nvtx,3);
v_aa = baryvec(:,:,aa,:);
s_bb = scalarC(:,:,bb,:);
v_cc = vectorC(:,:,cc,:);

if isequal(size(s_bb,[1,2]),[Ndim,Ndim])
	v_cc = permute(v_cc,[2,1,3,4]);
end
integrand = sum(v_aa .* s_bb .* v_cc,[1,2]);
integrand = permute(integrand,[3,4,1,2]);

% compute integration coefficients
svol = simplex_vol(pts,tri);
int_ikl = bary_integral(Ndim,4,2); % integration coeffs ignoring jj
int_jkl = bary_integral(Ndim,4,1); % integration coeffs ignoring ii

[~,ii,jj,kk,ll] = bary_integral(Ndim,4);
ii_jj = sub2ind([Nvtx,Nvtx],ii,jj);
ii_kk_ll = sub2ind([Nvtx,Nvtx,Nvtx],ii,kk,ll);
jj_kk_ll = sub2ind([Nvtx,Nvtx,Nvtx],jj,kk,ll);

% loop over simplexes
[rows,cols] = npermsk(Nvtx,2);
rows = tri(rows(:),:);
cols = tri(cols(:),:);
vals = zeros(size(rows));

for ss = 1:Ntri
	int_ss = integrand(:,ss);
	Xij = accumarray(ii_jj,int_ikl .* int_ss(jj_kk_ll,:) - int_jkl .* int_ss(ii_kk_ll,:));
	vals(:,ss) = svol(ss) .* Xij(:);
end
X = sparse(rows,cols,vals,Npts,Npts);
X = (X - X.')/2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Z] = femat_symm_vector(pts,tri,scalarC,vectorC)
% constructs the finite-element matrix for the real-valued vector-potential operator (2nd order)
arguments
	pts (:,:) double {mustBeFinite,mustBeReal}
	tri (:,:) {mustBeInteger,mustBePositive}
	scalarC double {mustBeFinite}
	vectorC double {mustBeFinite}
end
[Ndim,Npts] = size(pts);
[Nvtx,Ntri] = size(tri);

% construct simplex coefficients
scalarC = coeff_expand(pts,tri,scalarC);
vectorC = coeff_expand(pts,tri,vectorC);

% compute integrand values
[aa,bb,cc] = npermsk(Nvtx,3);
v_aa = vectorC(:,:,aa,:);
s_bb = scalarC(:,:,bb,:);
v_cc = vectorC(:,:,cc,:);

if isequal(size(s_bb,[1,2]),[Ndim,Ndim])
	v_cc = permute(v_cc,[2,1,3,4]);
end
integrand = sum(v_aa .* s_bb .* v_cc,[1,2]);
integrand = permute(integrand,[3,4,1,2]);

% compute integration coefficients
svol = simplex_vol(pts,tri);
[intCoeff,ii,jj,kk,ll,mm] = bary_integral(Ndim,5);
ii_jj = sub2ind([Nvtx,Nvtx],ii,jj);
kk_ll_mm = sub2ind([Nvtx,Nvtx,Nvtx],kk,ll,mm);

% loop over simplexes
[rows,cols] = npermsk(Nvtx,2);
rows = tri(rows(:),:);
cols = tri(cols(:),:);
vals = zeros(size(rows));

for ss = 1:Ntri
	int_ss = integrand(:,ss);
	Zij = accumarray(ii_jj, intCoeff .* int_ss(kk_ll_mm,:));
	vals(:,ss) = svol(ss) .* Zij(:);
end
Z = sparse(rows,cols,vals,Npts,Npts);
Z = (Z + Z.')/2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%