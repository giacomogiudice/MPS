function n = occupation_mps(occupation,state)
% Calculates occupation number for each site of a for given MPS of a density
% matrix. 
% Used for density-matrix MPS routines
%
% INPUT
%   N:  number of sites
%   d:  local Hilbert space dimension 
% OUTPUT
%   n:  array of corresponding occupation number for each site

N = length(state);
n= zeros(N,1);
for i = 1:N
    n(i) = braket(occupation{i},state);
end
end