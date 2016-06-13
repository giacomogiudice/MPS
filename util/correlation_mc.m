function g2 = correlation_mc(correlation,state,n)
% Calculates correlation matrix for given state. Each matrix element is
% computed as <a'a'aa>/(<a'a><a'a>)
% Used for 'monte_carlo'
%
% INPUT
%   N:  number of sites
%   d:  local Hilbert space dimension
%   n:  array of corresponding occupation number for each site
% OUTPUT
%   g2: correlation matrix (symmetric)

N = length(correlation);
g2 = zeros(N,N);
for i = 1:N
    for j = 1:i
        g2(i,j) = (state'*correlation{i,j}*state)/(real(n(i))*real(n(j)));
        g2(j,i) = g2(i,j);
    end
end
end