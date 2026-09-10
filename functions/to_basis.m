function [data] = to_basis(basis,data,shape,B)
% Converts the inputted data into the specified basis. 
% The inputted data can take the form of a matrix, or an array/pagearray containing
% row or column vectors. The input "B" corresponds to the basis conversion
% matrix for the desired transformation.
arguments
	basis string {mustBeMember(basis,["full","orth","sparse"])}
	data double {mustBeFinite}
	shape string {mustBeMember(shape,["row","col","mat"])}
	B (:,:) double {mustBeFinite}
end
ismat = (shape=="mat");
isrow = (shape=="row");
iscol = (shape=="col");

% ensure Q_rows > Q_cols
if size(B,1) < size(B,2); B=B.'; end
[Bmax, Bmin] = size(B);

% check size of input
if ismat
	if ~ismatrix(data) || size(data,1)~=size(data,2)
		error("ERROR: input is not a square matrix")
	elseif ~(all(size(data)==Bmax) || all(size(data)==Bmin))
		error("ERROR: the size of the input matrix does not agree with the basis-conversion matrix");
	end
elseif isrow && ~ismember(size(data,2),[Bmax,Bmin])
	error("ERROR: the size of the inputted row-vector does not agree with the basis-conversion matrix");
elseif iscol && ~ismember(size(data,1),[Bmax,Bmin])
	error("ERROR: the size of the inputted column-vector does not agree with the basis-conversion matrix");
end

% check if data is already in the correct basis
if basis=="full"
	if ismat && all(size(data)==Bmax); return; end
	if isrow && (size(data,2)==Bmax); return; end
	if iscol && (size(data,1)==Bmax); return; end
elseif ismember(basis,["orth","sparse"])
	if ismat && all(size(data)==Bmin); return; end
	if isrow && (size(data,2)==Bmin); return; end
	if iscol && (size(data,1)==Bmin); return; end
end

% determine if data is a pagearray
ispage = (ndims(data)==3 & all(size(data)~=0));

if ismat
	% check symmetry
	herm = ishermitian(data) - ishermitian(data,"skew");
	symm = issymmetric(data) - issymmetric(data,"skew");
end

% convert to basis
switch basis
	case "full"
		% convert
		if ismat
			data = ((B.') \ data / B);
		elseif iscol && ispage
			data = pagemtimes(B,data);
		elseif iscol & ~ispage
			data = B * data;
		elseif isrow & ispage
			data = pagemtimes(data,B);
		elseif isrow & ~ispage
			data = data * B;
		end

	case "orth"
		data = full(data);
		B = full(B);

		% convert
		if ismat
			data = (B.' * data * B);
		elseif iscol & ispage
			data = pagemtimes(inv(B),data);
		elseif iscol & ~ispage
			data = B \ data;
		elseif isrow & ispage
			data = pagemtimes(data,inv(B));
		elseif isrow & ~ispage
			data = data / B;
		end

	case "sparse"
		% convert
		if ismat
			data = (B.' * data * B);
		elseif iscol & ispage
			data = pagemtimes(B.',data);
		elseif iscol & ~ispage
			data = B.' * data;
		elseif isrow & ispage
			data = pagemtimes(data, B.');
		elseif isrow & ~ispage
			data = data * B.';
		end
end

% restore matrix symmetry
if ismat
	if herm==0 & symm==0
		% no symmetries
	elseif herm==-1
		% skew hermitian
		data = (data - data')/2;
	elseif herm==+1
		% hermitian
		data = (data + data')/2;
	elseif symm==-1
		% skew symmetric
		data = (data - data.')/2;
	elseif symm==+1
		% symmetric
		data = (data + data.')/2;
	end
end
% done
end