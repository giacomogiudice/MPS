function obs = expectationvalue(obs_set,mps)
% Computes the expectation values for a set of observables obs
% in MPO form.
%
% INPUT
%   obs_set:	rank-2 cell array corresponding to the set of observables
%				in MPO form. The first index corresponds to the ordering of
%				the observables, while the second index corresponds to the 
%				site (physical) index, and necessarily must match the number
%				of sites in mps
%   mps:      	cell array corresponding to the state to compute 
%				observables upon
% OUTPUT
%   obs:		vector corresponding to the observables associated to 
%				obs_set

[M,N] = size(obs_set);

obs = zeros(M,1);
for ind = 1:M
	block = 1;
	for site=N:-1:1
		block = update_right(block,mps{site},obs_set{ind,site});
	end
	obs(ind) = block;
end
end

function new_block = update_right(prev_block,M,O)
	new_block = contract(prev_block,1,M,2);
	new_block = contract(new_block,1,conj(M),2);
	new_block = contract(new_block,[2,4],O,[4,3]); 
end
