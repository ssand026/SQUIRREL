function [solverFunc] = setup_litsolve(A,dim,opt)
% Returns a function for iteratively solving a system of linear equations.
%
% Use if many linear solves are being done with the same large/sparse matrix.
% The method "gmres" works well for most matrices, or set method to "auto" to 
% automatically choose an iterative solver.
arguments
	A (:,:) double
	dim (1,1) double {mustBeMember(dim,[1,2])}
	opt.method string {mustBeMember(opt.method, ...
		["auto","lsqr","pcg","minres","gmres","qmr","tfqmr"])} = "gmres";
	opt.iterTol (1,1) double {mustBeInRange(opt.iterTol,0,1)} = 1e-9;
	opt.iterMax (1,1) double {mustBePositive,mustBeInteger} = 128;
end

% precondition using equilibrate
if dim==2; A = A.'; end
[P,R,C] = equilibrate(A);
A = R * (P * A) * C;

% pre-condition using an incomplete LU decomp
if issparse(A)
	[L,U] = ilu(A,struct('type','nofill','droptol',1e-1,'thresh',0));
else
	L = [];
	U = [];
end

if size(A,1)~=size(A,2)
	% A is rectangular
	opt.method = "lsqr"; 
elseif opt.method == "auto"
	% choose a solver method
	if issymmetric(A)
		try chol(A)
			opt.method = "pcg"; % A is positive definite
		catch
			opt.method = "minres"; % A is not positive definite
		end
	else
		opt.method = "gmres"; % A is square but non-symmetric
	end
end

% select solve function
tol = opt.iterTol;
maxit = opt.iterMax;
switch opt.method
	case "lsqr"
		itsolver = @(A,b) lsqr(A,b,tol,maxit,L,U);
	case "gmres"
		itsolver = @(A,b) gmres(A,b,[],tol,maxit,L,U);
	case "pcg"
		itsolver = @(A,b) pcg(A,b,tol,maxit,L,U);
	case "minres"
		itsolver = @(A,b) minres(A,b,tol,maxit,L,U);
	case "qmr"
		itsolver = @(A,b) qmr(A,b,tol,maxit,L,U);
	case "tfqmr"
		itsolver = @(A,b) tfqmr(A,b,tol,maxit,L,U);
end

solverFunc = @iterative_linsolve;
%----------------------------------------------------------
	function [x] = iterative_linsolve(b)
	if dim==2; b = b.'; end
	%
	b = (R*(P*b));
	x = zeros(size(b));
	for ii = 1:size(x,2)
		[temp,~] = itsolver(A,b(:,ii));
		x(:,ii) = temp;
	end
	x = C*x;
	%
	if dim==2; x = x.'; end
	end
%----------------------------------------------------------
end