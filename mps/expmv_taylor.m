function [y,flag,mv] = expmv_taylor(t,A,v,varargin)
% Compute exp(t*A)*v using a Taylor expansion
% The approximation of exp(t*A)*v is performed using only matrix-vector
% multiplications A*x. Hence, A does not have to be a matrix, but can be
% any linear operator, passed as a function handle. In this case, v does
% not have to strictly be a vector, but can be a tensor of any size. The
% output y will be a tensor of the same size.
%
% y = expmv_taylor(t,A,v)
% Perform the approximation using the default parameters. The time-step t
% can be real or complex, A is a matrix or an function handle, v is
% vector or a pseudo-vector that is applied to A.
%
% y = expmv_taylor(t,A,v,tol)
% Specify a halting tolerance on the residue at t only (t is assumed to be
% small).
%
% y = expmv_taylor(t,A,v,tol,m)
% Specify a tolerance and a maximum order m of the Taylor expansion. This
% parameter also corresponds to the maximum number of A*x multiplications.
%
% y = expmv_taylor(t,A,v,options)
% Specify options in a structure. The following options are supported, the
% default is specified in parentheses:
%   tol     - residue tolerance (1e-12)
%   m       - maximum number of orders to compute (100)
%   disp    - set to true to print out information at each iteration (0)
%
% [y,flag] = expmv_taylor(t,A,v,...)
% Also returns a convergence flag. If flag is 0, convergence was achieved,
% otherwise it is 1.
%
% [y,flag,mv] = expmv_taylor(t,A,v,...)
% Also returns the number mv of Taylor orders used.
%
% See also: LINALG/EXPMV.

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
m = min(m,n);       % Maximum order
normv = norm(v(:));

% Check if linear operator is a matrix of a function
if isnumeric(A)
    assert(size(A,2) == size(v,1),'Size mismatch between matrix and vector.');
    fapply = @(x) A*x;
elseif isa(A,'function_handle')
    fapply = @(x) A(x);
else
    error('Operator must be a matrix or a function handle.');
end

w = v/normv;
y = w;
for j = 1:m
    % Build the next step in the Taylor expansion
    w = (t/j)*fapply(w);
    y = y + w;
    resnorm = norm(w(:));
    if verbose
        fprintf('j = %d, resnorm = %.2e\n',j,resnorm);
    end
    % Stopping condition based on norm last term
    if resnorm <= tol
        break
    end
end
flag = 0;
mv = j;
if mv == m
    warning('EXMPV_TAYLOR:NoConvergence','No convergence within %d steps.',m);
    flag = 1;
end
y = normv*y;
end
