function n = occupation_mc(occupation,state)
% Calculates occupation number for each site for given state. 
% Used for 'monte_carlo'
%
% INPUT
%   N:  number of sites
%   d:  local Hilbert space dimension
% OUTPUT
%   n:  array of corresponding occupation number for each site

N = length(occupation);
n= zeros(N,1);
for i = 1:N
    n(i) = state'*occupation{i}*state;
end
end