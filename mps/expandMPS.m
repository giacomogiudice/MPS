function state = expandMPS(mps)
% expands an MPS into its equivalent representation as a quantum-mechanical
% state vector. This should only be used for debugging purposes on small
% systems.
% 
% INPUT
%	mps:	cell-array of rank-3 tensors representing the MPS. Indexing
%			convention is (bond,bond,physical)
% OUTPUT
%	state:	column-vector representation of the state

N = length(mps);
d = size(mps{1},3);

block = permute(mps{N},[1,3,2]);
for site = N-1:-1:1
	block = contract(block,N-site+1,N-site,mps{site},3,2);
end

state = reshape(squeeze(block),[d^N,1]);
end
