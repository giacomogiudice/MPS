function state = expand_mps(mps)
% expands an MPS into its equivalent representation as a quantum-mechanical
% state vector. This should only be used for debugging purposes on small
% systems.
% 
% INPUT
%   mps:    cell-array of rank-3 tensors representing the MPS. Indexing
%           convention is (bond,bond,physical)
% OUTPUT
%   state:  column-vector representation of the state

N = length(mps);
d = size(mps{1},3);
state = zeros(d^N,1);
% Compute the matrix of all possible combinations and store it in combs
c = cell(1,N);
[c{:}] = ndgrid(1:d);
combs = fliplr(cell2mat(cellfun(@(v)v(:),c,'UniformOutput',false)));
for pos = 1:d^N
    prod = 1;
    for i = 1:N
        prod = prod*mps{i}(:,:,combs(pos,i));
    end
    state(pos) = prod;
end

end

