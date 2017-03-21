function S = swapmatrix(d,N)
% returns the permutation matrix S that converts a 'reshape-first' density
% matrix to a 'combine-first' density matrix. Notice that the converse is
% given by S', since S is unitary
%
% INPUT
%   d:      local Hilbert space dimension
%   N:      number of sites
% OUTPUT
%   S:      permutation matrix


Z = sparse(d,d);
E = Z;
S = sparse(d^2^N,d^2^N);

c = cell(1,N);  
[c{:}] = ndgrid(1:d);
combs = fliplr(cell2mat(cellfun(@(v) v(:),c,'UniformOutput',false)));

col = 1;
for i = 1:length(combs)
    for j = 1:length(combs)
        S_prime = sparse(1);
        for site = 1:N
            E = Z;
            E(combs(j,site),combs(i,site)) = 1;
            S_prime = kron(S_prime,reshape(E,1,[]));
        end
        S(col,:) = S_prime;
        col = col + 1;
    end
end
end