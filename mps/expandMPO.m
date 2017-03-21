function operator = expandMPO(mpo)
% expands an MPO into its equivalent representation as a quantum-mechanical
% state vector. This should only be used for debugging purposes on small
% systems.
% 
% INPUT
%   mpo:        cell-array of rank-3 tensors representing the MPS. Indexing
%               convention is (bond,bond,physical,physical)
% OUTPUT
%   operator:  	matrix representation of the operator

N = length(mpo);
d = size(mpo{1},4);
operator = zeros(d^N,d^N);
% Compute the matrix of all possible combinations and store it in combs
c = cell(1,N);  
[c{:}] = ndgrid(1:d);
combs = fliplr(cell2mat(cellfun(@(v) v(:),c,'UniformOutput',false)));
for pos1 = 1:d^N
    for pos2 = 1:d^N
        prod = 1;
        for i = 1:N
            prod = prod*mpo{i}(:,:,combs(pos1,i),combs(pos2,i));
        end
        operator(pos1,pos2) = prod;
    end
end
end
