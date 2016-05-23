function mps_out = sweep(mps_in,mpo,direction,varargin)
% Computes a DMRG-style sweep on an MPS, applying some unitary operator
% and does canonization + decimation using 'canonize'
%
% INPUT
%   mps_in:         cell array corresponding to input MPS
%   mpo:            cell array corresponding to MPO to apply
%   direction:      specifies left (-1) or right (+1) canonization
%   D_max:          (optional) maximum bond size
%   epsilon:        (optional) maximum error in truncation
% OUTPUT
%   mps_out:        resulting MPS after computation
%   max_bond_flag:  returns true if the maximum bond length has been
%                   attained while compressing

N = length(mps_in);
mps_out = cell(1,N);

if isempty(varargin)
    D_max = [];
    epsilon = [];
else
    D_max = varargin{1};
    epsilon = varargin{2};
end
switch direction
    case +1 % Left sweep
        mps_out{1} = apply(mpo{1},mps_in{1});
        [mps_out{1},SV] = canonize(mps_out{1},1,D_max,epsilon);
        for site = 2:N
            mps_out{site} = apply(mpo{site},mps_in{site});
            mps_out{site} = contract(SV,2,mps_out{site},1);
            [mps_out{site},SV] = canonize(mps_out{site},1,D_max,epsilon);
        end
        mps_out{N} = mps_out{N}*sign(SV);
        
    case -1 % Right sweep
        mps_out{N} = apply(mpo{N},mps_in{N});
        [mps_out{N},US] = canonize(mps_out{N},-1,D_max,epsilon);
        for site = (N-1):(-1):1
            mps_out{site} = apply(mpo{site},mps_in{site});
            mps_out{site} = permute(contract(mps_out{site},2,US,1),[1 3 2]);
            [mps_out{site},US] = canonize(mps_out{site},-1,D_max,epsilon);
        end
        mps_out{1} = mps_out{1}*sign(US);
end
end
