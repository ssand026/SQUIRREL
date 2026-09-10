function [varargout] = npermsk(n,k)
% For a set containing "n" distinct elements, return all possible permutations of length "k"
% The order of the elements of a permutation does matter: each row of the output
% contains "k" elements in the range [1,n], where repeated values are allowed.
%
% SEE ALSO: NCOMBSK
arguments
	n (1,1) {mustBePositive,mustBeInteger}
	k (1,1) {mustBePositive,mustBeInteger}
end

% create the permutation list
perm_list = repmat([1:n^k].',1,k);
perm_list = ceil(perm_list ./ (n.^[0:k-1]));
perm_list = mod(perm_list,n);
perm_list(perm_list==0) = n;

if nargout<=1
	% output the list of permutations as an array
	varargout{1} = perm_list;
else
	% push each column of the permutation list to the output
	varargout = cell(1,nargout-1);
	for ii = 1:nargout
		varargout{ii} = perm_list(:,ii);
	end
end
end