function [mps_out,mps_norm] = sweep_fast(mps_in,direction)
% Computes a DMRG-style sweep on an MPS, canonizing it using QR
% decomposition in the direction specified.
%
% INPUT
%   mps_in:         cell array corresponding to input MPS
%   mpo:            cell array corresponding to MPO to apply
%   direction:      specifies left (-1) or right (+1) canonization
% OUTPUT
%   mps_out:        resulting MPS after computation
%   mps_norm:       square of the norm of the input MPS 

N = length(mps_in);
mps_out = cell(1,N);

switch direction
    case +1 % Left sweep
        [mps_out{1},R] = canonize_fast(mps_in{1},1);
        for site = 2:N
            mps_out{site} = contract(R,2,mps_in{site},1);
            [mps_out{site},R] = canonize_fast(mps_out{site},1);
        end
        
    case -1 % Right sweep
        [mps_out{N},R] = canonize_fast(mps_in{N},-1);
        for site = (N-1):(-1):1
            mps_out{site} = permute(contract(mps_in{site},2,R,1),[1 3 2]);
            [mps_out{site},R] = canonize_fast(mps_out{site},-1);
        end
end
mps_norm = R;
end
