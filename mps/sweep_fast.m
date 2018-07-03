function [mps_out,mps_norm] = sweep_fast(mps_in,direction)
% Computes a DMRG-style sweep on an MPS, canonizing it using QR
% decomposition in the direction specified.
%
% INPUT
%	mps_in:			cell array corresponding to input MPS
%	mpo:			cell array corresponding to MPO to apply
%	direction:		specifies left (-1) or right (+1) canonization
% OUTPUT
%	mps_out:		resulting MPS after computation
%	mps_norm:		norm of the output MPS 

N = length(mps_in);
mps_out = cell(1,N);

switch direction
	case +1 % Left sweep
		[mps_out{1},carryover] = canonize_fast(mps_in{1},1);
		for site = 2:N
			mps_out{site} = contract(carryover,2,2,mps_in{site},3,1);
			[mps_out{site},carryover] = canonize_fast(mps_out{site},1);
		end
		
	case -1 % Right sweep
		[mps_out{N},carryover] = canonize_fast(mps_in{N},-1);
		for site = (N-1):(-1):1
			mps_out{site} = permute(contract(mps_in{site},3,2,carryover,2,1),[1 3 2]);
			[mps_out{site},carryover] = canonize_fast(mps_out{site},-1);
		end
end
% Put the phase in the last term
mps_norm = abs(carryover);
mps_out{site} = sign(carryover)*mps_out{site};
end
