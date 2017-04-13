function sprod = braket(mps_1,mps_2)
% Returns the scalar product < MPS_1 | MPS_2 > that runs in polynomial time.
% The inputs are cell arrays corresponding to MPS decompositions,
% notice that the bond dimensions of the MPS does not have to be the same.
% 
% INPUT
%	mps_1, mps_2:	cell-array representations of MPSs to contract, each
%					element is a rank-3 tensor with index convention
%					(bond,bond,physical)
% OUTPUT
%	sprod:	scalar (hopefully, but not guaranteed) corresponding to
%			the scalar product

N = length(mps_1);
virtual_obs = cell(1,N);
sprod = sandwich(mps_1,virtual_obs,mps_2);
end
