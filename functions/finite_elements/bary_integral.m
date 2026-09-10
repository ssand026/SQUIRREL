function [intCoeffs,varargout] = bary_integral(Ndim, Nterm, ignore)
% Computes the barycentric integral coefficients
% 
% The output depends on the number of spatial dimensions and the total number of 
% integration terms. The optional argument "ignore" will exclude the corresponding 
% barycoordinate from the integral.
arguments
	Ndim  (1,1) double {mustBeInteger,mustBeNonnegative}
	Nterm (1,1) double {mustBeInteger,mustBePositive}
	ignore (1,:) double {mustBeInteger,mustBePositive} = [];
end

% construct the indices of the tensor
index = npermsk(Ndim+1,Nterm);
Nrows = size(index,1);

% push barycoord indices to the output
varargout = cell(1,nargout-1);
for jj = 1:nargout-1
	varargout{jj} = index(:,jj);
end

% remove any ignored terms from the index
inPoly = [1:Nterm];
if ~isempty(ignore)
	inPoly(any(inPoly==ignore(:),1)) = [];
end
index = index(:,inPoly);
polyDeg = length(inPoly);

% check if integral coefficients are uniform 
if polyDeg < 2
	intCoeffs = ones(Nrows,1) / factorial(Ndim + polyDeg); 
	return
end

% extract barycoord powers and compute coefficients
%----------------------------------------------------------
sorted = sort(index,2);       % sorted list of barycoord indices
gamma = cumprod([1:Nterm].'); % precompute the factorial values
exponents = ones(Nrows,1);
intCoeffs = ones(Nrows,1);

% loop over indices
for kk = 2:polyDeg
	% check entries where the referenced barycoord has changed
	isDiff = (sorted(:,kk-1)~=sorted(:,kk)); 

	% apply changed exponents to intCoeff, then update the exponents
	intCoeffs(isDiff) = intCoeffs(isDiff) .* gamma(exponents(isDiff)); 
	exponents = (~isDiff .* exponents) + 1;
end
intCoeffs = intCoeffs .* gamma(exponents);

% apply overall factor
intCoeffs = factorial(Ndim) .* intCoeffs / factorial(Ndim + polyDeg);
end

