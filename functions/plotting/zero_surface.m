function [zero_pts, zero_faces] = zero_surface(pts,tri,coeff)
% Interpolates over a mesh to find the zero-valued surface of a scalar field
arguments
	pts (:,:) double {mustBeFinite,mustBeReal}
	tri (:,:) {mustBeInteger,mustBePositive}
	coeff double {mustBeFinite,mustBeReal}
end
[Ndim,~] = size(pts);
tri = tri(1:Ndim+1,:);

% default to empty
zero_pts = [];
zero_faces = [];

% check if the coeff contains any zero-valued points
[minC,maxC] = bounds(coeff(:));
if (minC > 0) || (maxC < 0)
	return
end

% check/expand coefficient input
[coeff,order,cont_type] = coeff_check(pts,tri,coeff);
if order ~=0 
	error("ERROR: the input must be a scalar coefficient")
end
if ismember(cont_type,["domain-constant","simplex-constant"])
	return
else
	coeff = squeeze(coeff_expand(pts,tri,coeff));
end

% find any faces that lie exactly on the zero-valued surface
exact = logical(coeff == 0);
if nnz(exact) >= Ndim
	numExact = sum((coeff==0),1);

	facePerms = ncombsk(Ndim+1,Ndim);
	facePerms = reshape(facePerms.',[],1);
	
	% simplexes containing one zero-valued face
	oneFace = tri(:,numExact==Ndim);
	oneFace = reshape(oneFace(exact(:,numExact==Ndim)),Ndim,[]);
	
	% simplexes containing all zero-valued faces
	allFace = tri(:,numExact==Ndim+1);
	allFace = reshape(allFace(facePerms,:),Ndim,[]);

	zeroSurf = [oneFace, allFace];
	zeroSurf = reshape(pts(:,zeroSurf(:)),Ndim,Ndim,[]);
	zeroSurf = permute(zeroSurf,[2,3,1]);

	tri = tri(:,numExact<Ndim);
	coeff = coeff(:,numExact<Ndim);
else
	zeroSurf = [];
end

% check if the coeff contains any zero crossings
if (minC >= 0) || (maxC <= 0)
	% coeff does not cross zero
	if ~isempty(zeroSurf)
		[zero_pts,~,zero_faces] = unique(zeroSurf,"rows","stable");
		zero_faces = reshape(zero_faces,Ndim,[]).';
	end
	return
end

% get the subset of simplexes that span the zero surface
above = logical(coeff > 0);
below = logical(coeff < 0);

% check if the coeff contains any zero crossings
spansZero = any(above,1) & any(below,1);
if ~any(spansZero)
	% coeff does not cross zero
	if ~isempty(zeroSurf)
		[zero_pts,~,zero_faces] = unique(zeroSurf,"rows","stable");
		zero_faces = reshape(zero_faces,Ndim,[]).';
	end
	return
else
	Ntri = sum(spansZero);
	above = above(:,spansZero);
	below = below(:,spansZero);
	coeff = coeff(:,spansZero);
	tri = tri(:,spansZero);
end

% get all combinations of simplex edges and get the sets of orthogonal edges
[edg_a, edg_b] = ncombsk(Ndim+1,2);
Nedg = numel(edg_a & edg_b);

orthEdges = (edg_a'~=edg_b) & (edg_b'~=edg_a) & (edg_a'~=edg_a) & (edg_b'~=edg_b);
[row,col] = find(orthEdges);
orthPairs = unique(sort([row,col],2),"rows");

% for each set of simplex edges, find the zero point location, if it exists
zeroPts = NaN(Nedg,Ntri,Ndim);
zeroInd = false(Nedg,Ntri);
for ii = 1:Nedg
	doesCross = any(above([edg_a(ii),edg_b(ii)],:),1) & ...
				any(below([edg_a(ii),edg_b(ii)],:),1);
	
	edg_frac = coeff([edg_a(ii),edg_b(ii)],doesCross);
	edg_frac = abs(edg_frac) ./ abs(edg_frac(1,:)-edg_frac(2,:));

	zeroLoc = edg_frac(2,:) .* pts(:,tri(edg_a(ii),doesCross)) + ...
		      edg_frac(1,:) .* pts(:,tri(edg_b(ii),doesCross));

	zeroPts(ii,doesCross,:) = permute(zeroLoc,[3,2,1]);
	zeroInd(ii,:) = doesCross;
end
numPts = sum(zeroInd,1);

% get the set of zeroPts where the internal zero-surface forms a simplex
setA = zeroInd(:,(numPts==Ndim));
zeroPtsA = zeroPts(:,(numPts==Ndim),:);
surfA = zeroPtsA(repmat(setA,1,1,Ndim));
surfA = reshape(surfA,Ndim,[],Ndim);

% get the set of zeroPts where the internal zero-surface does NOT form a simplex
setB = zeroInd(:,(numPts==Ndim+1));
setC = zeroInd(:,(numPts==Ndim+1));
zeroPtsBC = zeroPts(:,(numPts==Ndim+1),:);

whereOrth = setB(orthPairs(:,1),:) & setC(orthPairs(:,2),:);
whereOrth = max([1:size(orthPairs,1)].' .* whereOrth,[],1);

b = orthPairs(whereOrth,1).';
setB(sub2ind(size(setB),b,1:length(b))) = false;
surfB = zeroPtsBC(repmat(setB,1,1,Ndim));
surfB = reshape(surfB,Ndim,[],Ndim);

c = orthPairs(whereOrth,2).';
setC(sub2ind(size(setC),c,1:length(c))) = false;
surfC = zeroPtsBC(repmat(setC,1,1,Ndim));
surfC = reshape(surfC,Ndim,[],Ndim);

% combine all interpolated zero-surfaces into a unified set
zeroSurf = [zeroSurf, surfA, surfB, surfC];
zeroSurf = reshape(zeroSurf,[],Ndim);

[zero_pts,~,zero_faces] = unique(zeroSurf,"rows","stable");
zero_faces = reshape(zero_faces,Ndim,[]).';
end