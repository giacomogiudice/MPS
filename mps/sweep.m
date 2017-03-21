function [mps_out,mps_norm] = sweep(mps_in,mpo,direction,varargin)
% Computes a DMRG-style sweep on an MPS, applying some operator in MPO form
% and then doing canonization + decimation using 'canonize'
%
% INPUT
%   mps_in:         cell array corresponding to input MPS
%   mpo:            cell array corresponding to MPO to apply
%   direction:      specifies left (-1) or right (+1) canonization
%   D_max:          (optional) maximum bond size
%   epsilon:        (optional) maximum error in truncation
% OUTPUT
%   mps_out:        resulting MPS, in the opposite canonization of the
%                   input MPS
%   mps_norm:       square of the norm of the input MPS

if isempty(varargin)
    if isempty(mpo)
        [mps_out,mps_norm] = sweep_fast(mps_in,direction);
        return
    end
    D_max = [];
    epsilon = [];
else
    D_max = varargin{1};
    epsilon = varargin{2};
end

N = length(mps_in);

if isempty(mpo)
    mult = @(W,M) M;
    mpo = cell(1,N);    
else
    mult = @apply;
end

mps_out = cell(1,N);

switch direction
    case +1 % Left sweep
        mps_out{1} = mult(mpo{1},mps_in{1});
        [mps_out{1},SV] = canonize(mps_out{1},1,D_max,epsilon);
        for site = 2:N
            mps_out{site} = mult(mpo{site},mps_in{site});
            mps_out{site} = contract(SV,2,mps_out{site},1);
            [mps_out{site},SV] = canonize(mps_out{site},1,D_max,epsilon);
        end
        mps_norm = SV;
        
    case -1 % Right sweep
        mps_out{N} = mult(mpo{N},mps_in{N});
        [mps_out{N},US] = canonize(mps_out{N},-1,D_max,epsilon);
        for site = (N-1):(-1):1
            mps_out{site} = mult(mpo{site},mps_in{site});
            mps_out{site} = permute(contract(mps_out{site},2,US,1),[1 3 2]);
            [mps_out{site},US] = canonize(mps_out{site},-1,D_max,epsilon);
        end
        mps_norm = US;
    end
end
end
