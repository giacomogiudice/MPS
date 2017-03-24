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
	obs(ind) = sandwich(mps,obs_set(ind,:),mps);
end
end
