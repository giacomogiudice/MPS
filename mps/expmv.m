function [y,flag,mv] = expmv(t,A,v,varargin)
% Compute exp(t*A)*v using a Krylov subspace method.
% The approximation of exp(t*A)*v is performed using only matrix-vector
% multiplications A*x. Hence, A does not have to be a matrix, but can be
% any linear operator, passed as a function handle. In this case, v does
% not have to strictly be a vector, but can be a tensor of any size. The
% output y will be a tensor of the same size.
%
% y = expmv(t,A,v)
% Perform the approximation using the default parameters. The time-step t
% can be real or complex, A is a matrix or an function handle, v is
% vector or a pseudo-vector that is applied to A.
%
% y = expmv(t,A,v,tol)
% Specify a halting tolerance on the residue at t only (t is assumed to be
% small).
%
% y = expmv(t,A,v,tol,m)
% Specify a tolerance and a maximum size m of the Krylov subspace. This
% parameter also corresponds to the maximum number of A*x multiplications.
%
% y = expmv(t,A,v,options)
% Specify options in a structure. The following options are supported, the
% default is specified in parentheses:
%   tol     - residue tolerance (1e-12)
%   m       - maximum number of Krylov vectors to store (100)
%   disp    - set to true to print out information at each iteration (0)
%
% [y,flag] = expmv(t,A,v,...)
% Also returns a convergence flag. If flag is 0, convergence was achieved,
% otherwise it is 1.
%
% [y,flag,mv] = expmv(t,A,v,...)
% Also returns the number mv of Krylov vector used.
%
% This is an adapted version of the algorithm presented in the following
% references:
%   R. B. Sidje, Expokit: A software package for computing matrix
% exponentials (1998)
%   M. A. Botchev, A short guide to exponential Krylov subspace time
% integration for Maxwell's equations (2012)

% Default parameters
tol = 1e-12;
m = 100;
verbose = false;

% Handle optional arguments
switch length(varargin)
    case 1
        if isstruct(varargin{1})
            options = varargin{1};
            if isfield(options,'tol')
                tol = options.tol;
            end
            if isfield(options,'m')
                m = options.m;
            end
            if isfield(options,'disp')
                verbose = options.disp;
            end
        else
            tol = varargin{1};
        end
    case 2
        tol = varargin{1};
        m = varargin{2};
end

s = size(v);        % Size of input tensor
n = prod(s);        % Number of elements of input tensor
m = min(m,n);       % Krylov space size
V = zeros(n,m+1);   % Projector onto Krylov space
H = zeros(m+1,m);   % Hessenberg matrix
normv = norm(v(:));
V(:,1) = v(:)/normv;

% Check if linear operator is a matrix of a function
if isnumeric(A)
    assert(size(A,2) == size(v,1),'Size mismatch between matrix and vector.');
    fapply = @(x) A*x;
elseif isa(A,'function_handle')
    fapply = @(x) reshape(A(reshape(x,s)),[n,1]);
else
    error('Operator must be a matrix or a function handle.');
end

for j = 1:m
    % Build Krylov subspace and Gram-Schmidt orthogonalize
    w = fapply(V(:,j));
    for i = 1:j
        H(i,j) = w'*V(:,i);
        w = w - H(i,j)*V(:,i);
    end
    H(j+1,j) = norm(w);
    % Build exponential approximation in Krylov subspace
    U = expm(t*H(1:j,1:j));
    resnorm = norm(-H(j+1,j)*U(j,1));
    if verbose
        fprintf('j = %d, resnorm = %.2e\n',j,resnorm);
    end
    % Stopping condition based on norm of residual
    if resnorm <= tol
        break
    end
    V(:,j+1) = w/H(j+1,j);
end
flag = 0;
mv = j;
if mv == m && resnorm > tol
    warning('EXMPV:NoConvergence','No convergence within %d steps.',m);
    flag = 1;
end
% Build resulting approximation in full space
y = normv*reshape(V(:,1:j)*U(:,1),s);
end
