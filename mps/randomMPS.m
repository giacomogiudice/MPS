function mps = randomMPS(N,D_max,d,direction,isreal)
% Generates a random mps of length N, bond length D_max and local H-space d
% then canonizes according to direction specified: left (+1) or right (-1).
% The maximal bond dimension at each site assumes open boundary conditions.
%
% INPUT
%	N:			number of sites
%	D_max:		maximum bond dimension of each size
%	d:			local Hilbert space dimension
%	direction:	(optional) specifies left (-1) or right (+1) canonization. 
%				Any other value implies no canonization
%	isreal		(optional) boolean to set to true for real-valued MPS.
%				Complex-valued by default
% OUTPUT
%	mps:	resulting random MPS 

if nargin < 5
	isreal = false;
end
mps = cell(1,N);
% Build matrix containing size of bond dimensions for the MPS
n = 0:(N-1);
bnd = min([min(d.^n,d.^(N-n));min(d.^(n+1),d.^(N-n-1))],D_max);

for site = 1:N
	mps{site} = randn([bnd(1,site),bnd(2,site),d]) + ~isreal*1i*randn([bnd(1,site),bnd(2,site),d]);
end

if nargin >= 4 && (direction == 1 || direction == -1)
	mps = sweep(mps,{},direction);
end
end
