function a = annihilation(N,d)
% Initializes annihilation operators in d^N Hilbert space
%
% INPUT
%   N:  number of sites
%   d:  local Hilbert space dimension
% OUTPUT
%   a:  cell array of annihilation operators, one for each site

a = cell(1,N);
A = sparse(diag(sqrt(1:(d-1)),1));
I = speye(d);
for k = 1:N
    a{k} = sparse(1);
    for j = 1:N
        if j == k
            a{k} = kron(a{k},A);
        else
            a{k} = kron(a{k},I);
        end
    end
end