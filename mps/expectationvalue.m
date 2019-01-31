function obs = expectationvalue(mps,obs_set,bounds)
% Computes expectation values of local observables for an MPS.
% The behavior depends on which observables is desired. 

% INPUT
% 	mps:		cell array correponding to MPS representation of the state
% 				WARNING: Assumes right (-1) canonization
%	obs_set:	behavior depends on what is provided. Possible options are:
%				'correlationlength':	correlation length computed as the
%										second largest eigenvalue of the 
%										transfer matrix
%				{M}:					expectation value of a one-point
%										observable M (d x d matrix)
%				{M_1,M_2}:				correlation <M_1 M_2> in the 
%										interval specified by the bounds
%										each operator is a (d x d) matrix
%	bounds:		(optional) array corresponding to the initial and the final
%				site of the desired interval
% OUTPUT
%	obs:		array corresponding to the measured observables

N = length(mps);
if nargin <= 2
	N_start = 1;
	N_end = N;
else
	N_start = bounds(1);
	N_end = bounds(2);
end
assert((N_start < N_end) & (N_end <= N),'Bounds provided are not valid.');
ind = 1;

if isnumeric(obs_set)
	obs_set = {obs_set};
end

if ischar(obs_set)	
	if strcmp(obs_set,'correlationlength')
		obs = zeros(1,N_end-N_start+1);
		options.issym = false;
		options.isreal = false;
		for site = N_start:N_end
			[D1,D2,d] = size(mps{site});
			if D1 ~= D2
				obs(ind) = NaN;
				ind = ind + 1;
				continue
			end
			D = D1;
			id = reshape(eye(D),[],1)/sqrt(D);
			applyT = @(v) reshape(update_block(reshape(v,[D,D]),mps{site},{},mps{site},-1),[D^2,1]) - (id'*v)*id;
			lambda = eigs(applyT,D^2,1,'lm',options);
			obs(ind) = -1/log(abs(lambda));
			ind = ind + 1;
		end
	else
		error('Unrecognized observable type %s.',obs_set);
	end
elseif iscell(obs_set)
	% Get fixed point up to N_start
	rho = 1;
	for site = 1:(N_start-1)
		rho = update_block(rho,mps{site},{},mps{site},1);
	end
	rho = rho/trace(rho);
	if length(obs_set) == 1
		% One-point observables
		obs = zeros(1,N_end-N_start+1);
		op = reshape(obs_set{1},[1,1,size(obs_set{1})]);
		for site = N_start:N_end
			rho_op = update_block(rho,mps{site},op,mps{site},1);
			obs(ind) = trace(rho_op);
			rho = update_block(rho,mps{site},{},mps{site},1);	
			ind = ind + 1;
		end
	elseif length(obs_set) == 2
		% Two-point observables
		obs = zeros(1,N_end-N_start);
		op1 = reshape(obs_set{1},[1,1,size(obs_set{1})]);
		op2 = reshape(obs_set{2},[1,1,size(obs_set{2})]);
		rho_op = update_block(rho,mps{N_start},op1,mps{N_start},1);
		for site = (N_start+1):N_end
			obs(ind) = trace(update_block(rho_op,mps{site},op2,mps{site},1));
			rho_op = update_block(rho_op,mps{site},{},mps{site},1);
			ind = ind + 1;
		end
	else
		error('Only one-point and two-point correlators are supported.');
	end
end
end
