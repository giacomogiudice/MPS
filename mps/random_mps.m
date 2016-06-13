function mps = random_mps(N,D_b,d,direction)
% Generates a random mps of length N, bond length D_b and local H-space d
% then canonizes according to direction specified: left (+1) or right (-1)
%
% INPUT
%   N:              number of sites
%   D_b:            bond dimension of each size, notice that canonization 
%                   may change this 
%   d:              local Hilbert space dimension
%   direction:      specifies left (-1) or right (+1) canonization. Any
%                   other value will imply no canonization
% OUTPUT
%   mps:            resulting random MPS 
%              
mps = cell(1,N);
mps{1} = rand(1,D_b,d);
for i = 2:N-1
    mps{i} = rand(D_b,D_b,d);
end
mps{N} = rand(D_b,1,d);

switch direction
    case +1 % Left sweep
        [mps{1},SV] = canonize(mps{1},1);
        for site = 2:N
            mps{site} = contract(SV,2,mps{site},1);
            [mps{site},SV] = canonize(mps{site},1);
        end
        mps{N} = mps{N}*sign(SV);
    case -1 % Right sweep
        [mps{N},US] = canonize(mps{N},-1);
        for site = (N-1):(-1):1
            mps{site} = permute(contract(mps{site},2,US,1),[1 3 2]);
            [mps{site},US] = canonize(mps{site},-1);
        end
        mps{1} = mps{1}*sign(US);
end
end