function [pts,tri] = mesh_anneal(pts,tri,opt)
% Improves the quality of the mesh by selectively swapping connected edges.
% Only works for 2D meshes.
arguments
	pts (2,:) double {mustBeFinite,mustBeReal}
	tri (:,:) {mustBeInteger,mustBePositive}
	opt.symmIndex (:,1) {mustBeInteger,mustBePositive} = [];
	opt.maxIters (1,1) {mustBeInteger,mustBePositive} = ceil(size(tri,2)/9);
end

if isempty(opt.symmIndex)
	S = 1; 
else 
	S = length(unique(opt.symmIndex));
end

% find the ideal number of triangles for each element
[a1,a2,a3] = mesh_angles(pts,tri);
total_ang = full(sparse(reshape(tri(1:3,:).',[],1),1,[a1; a2; a3],size(pts,2),1));
ideal_tri = round(total_ang/(pi/3),3);
clear a1 a2 a3 total_ang

% optimization loop
%-------------------------------------------------
quad_hist = [0 0 0 0];
iter = 1;
while iter < opt.maxIters
	% get mesh connectivity
	[pts,tri] = mesh_smooth(pts,tri,7);
	[pt_conn,num_tri,bnd_indx,~] = mesh_conn(pts,tri);

	% find adjacent nodes and triangle indices for quad regions
	pt_conn(bnd_indx) = 0;
	[row,col,val] = find(pt_conn);

	% ensure adjacent points have correct dimensions
	ind = (tri(1:3,val)~=([row,row,row].') & tri(1:3,val)~=([col,col,col].'));
	row = row(sum(ind,1)==1);
	col = col(sum(ind,1)==1);
	val = val(sum(ind,1)==1);
	adj = tri(1:3,val);
	adj = adj(adj~=([row,row,row].') & adj~=([col,col,col].'));
	
	% info about the edges pts(:,row)-pts(:,col) and the connected triangles
	quads = sortrows([sort([row,col],2),adj,val]);
	if isempty(quads); return; end

	% exclude any lone edges
	test = diff(quads(:,[1 2]),1,1);
	test = all(~test,2);
	test = (([test; 0] + [0; test])==1);
	quads = quads(test,:);
	if isempty(quads); return; end

	% restructure quad-list
	quads = [quads(1:2:end,[1 2 3 4]), quads(2:2:end,[3 4])];
	quads = quads(:,[1 2 3 5 4 6]);

	% ensure quads fall within a single geometric region
	same_region = reshape(tri(4,quads(:,5))==tri(4,quads(:,6)),[],1);
	if ~isempty(opt.symmIndex)
		symReg1 = opt.symmIndex(quads(:,5)); symReg1 = symReg1(:);
		symReg2 = opt.symmIndex(quads(:,6)); symReg2 = symReg2(:);

		same_region = same_region & (symReg1==symReg2) & (imag(symReg1)==0) & (imag(symReg2)==0);
		quads = [quads, symReg1, symReg2];
	end
	quads = quads(same_region,:);
	if isempty(quads); return; end

	% find locations where removing quads is an improvement to connectivity
	ideal = reshape(ideal_tri(quads(:,[1 2 3 4])),[],4);
	ntri1 = reshape(num_tri(quads(:,[1 2 3 4])),[],4);
	ntri2 = ntri1 + [-1, -1, 1, 1];
	Etri1 = ((ntri1-ideal).^2);
	Etri2 = ((ntri2-ideal).^2);

	% select for quadrilaterals where swapping the internal edge decreases mesh energy
	% edge lengths are used as a symmetry identifier (might need to swap to explicit indexing)
	Ediff = (Etri2 - Etri1) * [1;1;1;1];
	edge_len = -sqrt(sum((pts(:,quads(:,1)) - pts(:,quads(:,2))).^2,1).');

	quads = sortrows([Ediff,edge_len,quads]);
	quads = quads((quads(:,1)<0),:);
	if isempty(quads); return; end
	
	if size(quads,1) < S
		swap = [];
	else
		% find first set of quads with consistent energy/edge_len
		tol = 3*max(max(eps(pts(:,quads(:,[3 4])))));
		for qi = 0:(size(quads,1)-S)
			swap = qi + [1:S];
			dist = quads(swap,2);
			if all(dist>=(dist(1)-tol) & dist<=(dist(1)+tol)); break; end
			if qi==(size(quads,1)-S); return; end
		end
	end
	if max(swap) > size(quads,1); return; end

	% re-configure quads that improve energy
	quads(:,[1 2])=[];
	tri(1:3,quads(swap,5)) = quads(swap,[1 3 4]).';
	tri(1:3,quads(swap,6)) = quads(swap,[2 3 4]).';

	% check if oscillating
	quad_hist = [quad_hist(2:4),size(quads,1)];
	if quad_hist(1)==quad_hist(3) && quad_hist(2)==quad_hist(4) && quad_hist(2)~=quad_hist(3)
		return;
	end

	% update iter number
	iter = iter + 1;
end
end