function mps = randomMPS(N,D_b,d,direction)
% Generates a random mps of length N, bond length D_b and local H-space d
% then canonizes according to direction specified: left (+1) or right (-1).
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
mps{1} = randn(1,D_b,d) + 1i*randn(1,D_b,d);
for i = 2:N-1
    mps{i} = randn(D_b,D_b,d) + 1i*randn(D_b,D_b,d);
end
mps{N} = randn(D_b,1,d) + 1i*randn(D_b,1,d);

if direction == 1 || direction == -1
	mps = sweep(mps,{},direction);
end
end