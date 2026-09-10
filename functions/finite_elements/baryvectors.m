function [bvec] = baryvectors(pts,tri)
% (memoized) Computes the barycentric vector corresponding to each vertex in the mesh
arguments
	pts (:,:) double {mustBeFinite,mustBeReal}
	tri (:,:) {mustBeInteger,mustBePositive}
end

% enable memoization for repeated calls
persistent memFunc
if isempty(memFunc)
	memFunc = memoize(@baryvectors_main);
	memFunc.CacheSize = 1;
end
bvec = memFunc(pts,tri);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [bvec] = baryvectors_main(pts,tri)
% computes the set of barycentric vectors for each simplex in the mesh
[Ndim,~] = size(pts);
[~,Ntri] = size(tri);

% rearrange mesh data
mesh = zeros(Ndim,Ntri,Ndim+1);
for jj = 1:Ndim+1
	mesh(:,:,jj) = pts(:,tri(jj,:)');
end
mesh = permute(mesh,[1,3,2]);

% construct baryvectors
bvec = zeros(size(mesh));
for kk = 1:(Ndim+1)
	ind = [1:kk-1,kk+1:Ndim+1];

	% get side-length vectors
	svec = mesh(:,ind,:) - mesh(:,ind(end),:);
	svec(:,end,:) = 1;
	nS = vecnorm(svec,2,1);
	svec = svec ./ (nS + ~nS);

	% get normal vector (via Gram-Schmidt orthogonalization)
	Q = zeros(size(svec));
	for jj = 1:size(Q,2)
		Q(:,jj,:) = svec(:,jj,:);
		if jj > 1
			proj = sum(Q(:,1:jj-1,:).*Q(:,jj,:),1);
			proj = sum(Q(:,1:jj-1,:).*proj,2);
			Q(:,jj,:) = Q(:,jj,:) - proj;
		end
		% normalize
		nQ = vecnorm(Q(:,jj,:),2,1);
		Q(:,jj,:) = Q(:,jj,:) ./ (nQ + ~nQ);
	end
	nvec = Q(:,end,:);

	% compute the barycentric vectors
	lvec = mesh(:,kk,:) - mesh(:,ind(1),:);
	bvec(:,kk,:) = nvec ./ sum(nvec .* lvec,1);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%