function [out] = expval(bra,op,ket)
% Compute the expectation values for the arrays/pagearrays of bra/kets
arguments
	bra (:,:,:) double {mustBeFinite}
	op  (:,:)   double {mustBeFinite}
	ket (:,:,:) double {mustBeFinite}
end
Nb = ndims(bra);
Nk = ndims(ket);
Nd = max(size(bra,3),size(ket,3));

if ~isscalar(op)
	if isempty(op)
		% no value was given for the operator
		op = 1;
	elseif size(op,1)~=size(bra,2) | size(op,2)~=size(ket,1)
		% the operator does not commute with the bra/ket
		error("ERROR: sizes of the operator and the bra/ket do not commute")
	elseif isdiag(op) && all(abs(diag(op)-op(1,1)) <= 2*eps(op(1,1)))
		% the operator is an identity matrix
		op = op(1,1);
	end
end

% initialize output
out = zeros(size(bra,1),size(ket,2),Nd);

% check if operator is zero-valued
if all(op==0)
	return
end

% compute the expectation values
if Nb==2 & Nk==2 
	% bra and ket are matrices
	out = bra * op * ket;
elseif (Nb==2 & Nk==3)
	% ket is a pagearray
	bra = bra * op;
	for jj = 1:size(ket,3)
		out(:,:,jj) = bra * ket(:,:,jj);
	end
elseif (Nb==3 & Nk==2) 
	% bra is a pagearray
	ket = op * ket;
	for jj = 1:size(bra,3)
		out(:,:,jj) = bra(:,:,jj) * ket;
	end
elseif (Nb==3 & Nk==3) 
	% bra and ket are pagearrays
	if size(bra,3)~=Nd | size(ket,3)~=Nd
		error("ERROR: bra and ket have a different number of pages")
	end
	for jj = 1:Nd
		out(:,:,jj) = bra(:,:,jj) * op * ket(:,:,jj);
	end
else
	% bra and kets have invalid sizes
	error("ERROR: expval not available for given inputs")
end
end