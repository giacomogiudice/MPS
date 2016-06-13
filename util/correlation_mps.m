function g2 = correlation_mps(correlation,state,n)
% Calculates correlation matrix for a given MPS encoding a density matrix. 
% Each matrix element is computed as <a'a'aa>/(<a'a><a'a>)
% Used for density-matrix MPS routines
%
% INPUT
%   N:  number of sites
%   d:  local Hilbert space dimension
%   n:  array of corresponding occupation number for each site
% OUTPUT
%   g2: correlation matrix (symmetric)

N = length(state);
g2 = zeros(N,N);
for i = 1:N
    for j = 1:i
        g2(i,j) = braket(correlation{i,j},state)/(real(n(i))*real(n(j)));
        g2(j,i) = g2(i,j);
    end
end
end