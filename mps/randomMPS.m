function mps = randomMPS(N,D_max,d,direction)
% Generates a random mps of length N, bond length D_max and local H-space d
% then canonizes according to direction specified: left (+1) or right (-1).
% The maximal bond dimension at each site assumes open boundary conditions.
%
% INPUT
%	N:			number of sites
%	D_max:		maximum bond dimension of each size
%	d:			local Hilbert space dimension
%	direction:	specifies left (-1) or right (+1) canonization. Any other
%				value will imply no canonization
% OUTPUT
%	mps:	resulting random MPS 

mps = cell(1,N);
% Build matrix containing size of bond dimensions for the MPS
bnd = zeros(2,N);
first_half = floor(N/2);
bnd(:,1:first_half) = [d.^((1:first_half)-1);d.^(1:first_half)];
if mod(N,2)
	first_half = first_half + 1;
	bnd(:,first_half) = [d^(first_half-1);d^(first_half-1)];	
end
first_half = first_half + 1;
bnd(:,first_half:N) = [d.^(N+1-(first_half:N));d.^(N-(first_half:N))];
bnd = min(bnd,D_max);

for site = 1:N
	mps{site} = randn(bnd(1,site),bnd(2,site),d) + 1i*randn(bnd(1,site),bnd(2,site),d);
end

if direction == 1 || direction == -1
	mps = sweep(mps,{},direction);
end
end
