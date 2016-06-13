function psi = ground(N,d)
% Initializes a ground state in d^N Hilbert space
%
% INPUT
%   N:      number of sites
%   d:      local Hilbert space dimension
% OUTPUT
%   psi:    ground state

vac = zeros(d, 1);
vac(1) = 1;
psi = 1;
for k = 1:N
    psi = kron(psi,vac);
end
end
