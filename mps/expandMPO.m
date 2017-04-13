function operator = expandMPO(mpo)
% expands an MPO into its equivalent representation as a quantum-mechanical
% state vector. This should only be used for debugging purposes on small
% systems.
% 
% INPUT
%	mpo:		cell-array of rank-3 tensors representing the MPS. Indexing
%				convention is (bond,bond,physical,physical)
% OUTPUT
%	operator:	matrix representation of the operator

N = length(mpo);
d = size(mpo{1},3);
block = permute(mpo{N},[1,3,4,2]);
for site = N-1:-1:1
	block = contract(block,2*(N-site)+1,2*(N-site)-1,mpo{site},4,2);
end

operator = reshape(permute(squeeze(block),[1:2:2*N,2:2:2*N]),[d^N,d^N]);
end
