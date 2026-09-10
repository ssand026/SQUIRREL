function [solverFunc] = setup_linsolve(A,dim,method)
% Returns a function for solving the system of linear equations A*x = b.
% This prevents repeated calculation of any preconditioning steps. For small
% systems where extreme accuracy is not needed, use the method "preinvert";
% otherwise, use the method "decompose".
arguments
	A (:,:) double {mustBeFinite}
	dim (1,1) double {mustBeMember(dim,[1,2])} = 1;
	method string {mustBeMember(method,["decompose","preinvert","full"])} = "decompose"
end
switch method
	case "decompose"
		dA = decomposition(A);
		if dim==1; solverFunc = @(b) dA \ b; end
		if dim==2; solverFunc = @(b) b / dA; end
	case "preinvert"
		iA = inv(A);
		if dim==1; solverFunc = @(b) iA * b; end
		if dim==2; solverFunc = @(b) b * iA; end
	case "full"
		% not recommended
		if dim==1; solverFunc = @(b) A \ b; end
		if dim==2; solverFunc = @(b) b / A; end
end
end