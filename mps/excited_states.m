function [mps_list,energies,iters,variances] = excited_states(guesses,mpo,M,varargin)
% This function computes the lowest M eigenstates of a given MPO.
% The starting point for each minimization is initialized at random.
%
% INPUT
%	guesses:		cell array of initial MPS to provide to the optimizerA
%	mpo:			cell array corresponding to MPO, assumed to be Hermitian
%	M:				number of eigenstates to compute     
%	varargin:		(optional) options to pass to 'projected_search'
% OUTPUT
%	mps_list:		cell array of resulting MPS eigenstates
%	energies:		vector of eigenvalues of corresponding eigenstates		
%	iters:			number of iterations used for each optimization
assert(M >= 1,'Number of desired states must be larger than zero.');
N = length(mpo);
d = size(mpo{1},4);

constr_list = {};
mps_list = {};
energies = zeros(M,1);
iters = zeros(M,1);
variances = zeros(M,1);

for m = 1:M
	fprintf('Computing state #%d...\n',m);
	[state,E,iter] = projected_search(guesses{m},mpo,constr_list,varargin{:});
	mps_list{m} = state;
	constr_list{m} = sweep(state,{},1);
	energies(m) = E;
	iters(m) = iter;
	variances(m) = error_variance(state,mpo);
end
fprintf('Done\n');
end