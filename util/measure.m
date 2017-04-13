function obs = measure(obs_set,density_mps,obs_type)
% Calculates the values associated to each observable for a density matrix.
% Both the density matrix and each observable have to be in MPDO form, with
% the two physical indices bunched together in the last index.
% NOTE: result should be rescaled depending on the normalization
% of the density matrix.
%
% INPUT
%	obs_set:		rank-1 or rank-2 cell array, where each element is 
%					an MPS corresponding to the observable
%	density_mps:	MPS corrensponding to the MPDO
%	obs_type:		string corresponding to the type of observable:
%					'single' for N sigle sites, 'neighbor' for N-1 two-site
%					nearest-neighbor observables, and 'set' to compute
%					a 2D set.
% OUTPUT
%	obs:  			array of corresponding observable for each site

M = length(density_mps);

switch obs_type
case 'single'
	obs = zeros(N,1);
	for site = 1:N
		obs(site) = braket(obs_set{site},density_mps);
	end
case 'neighbor'
	obs = zeros(N-1,1);
	for site = 1:(N-1)
		obs(site) = braket(obs_set{site},density_mps);
	end
case 'set'
	obs = zeros(M,N);
for i = 1:size(obs_set,1)
	for j = 1:size(obs_set,2)
		obs(i,j) = braket(obs_set{i,j},density_mps);
	end
end
otherwise
	fprintf('Unrecognized observable type: %s\n',obs_type);
	obs = [];
end
end
