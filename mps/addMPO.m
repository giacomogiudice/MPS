function mpo_out = addMPS(mpo_1,mpo_2)
% Adds two MPO by combining each element in a block-diagonal way
%
% INPUT
%   mpo_1,mpo_2:    MPO represented as a cell array, with
%                   index convention (bond,bond,physical,physical)
% OUTPUT
%   mpo_2:          MPO corresponding to the sum of the inputs

N = length(mpo_1);
assert(length(mpo_2) == N,'Mismatch between sizes of input MPOs.')

mps_1 = cellfun(@(c) reshape(c,[size(c,1),size(c,2),size(c,3)*size(c,4)]),mpo_1,'UniformOutput',false);
mps_2 = cellfun(@(c) reshape(c,[size(c,1),size(c,2),size(c,3)*size(c,4)]),mpo_2,'UniformOutput',false);
mps_out = addMPS(mps_1,mps_2);

for site = 1:N
    [b_1,b_2,~] = size(mps_out{site});
    [~,~,d_1,d_2] = size(mpo_1{site});
    mpo_out{site} = reshape(mps_out{site},[b_1,b_2,d_1,d_2]);
end
end
