function block = sandwich(mps_1,mpo,mps_2)
% Returns the matrix element < MPS_1 | MPO | MPS_2 > that runs in
% polynomial time.
%
% INPUT
%   mps_1, mps_2:   cell array representations of MPSs to contract, each
%                   element is a rank-3 tensor with index convention
%                   (bond,bond,physical)
%	mpo:			cell array corresponding to the operator to contract in
%					between the states
% OUTPUT
%   block:          scalar corresponding to the matrix element

N = length(mps_1);
block = 1;

for site=N:-1:1
	block = update_block(block,mps_1{site},mpo{site},mps_2{site},-1);
end
end
