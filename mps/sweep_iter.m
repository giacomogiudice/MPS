function [mps_out,mps_norm,iter] = sweep_iter(mps_in,mpo,mps_out,iter_max,tolerance,max_storage_size)
% Computes a DMRG-style sweep on an MPS, applying some operator in MPO form
% and then doing canonization and decimation using the iterative method.
% WARNING: The MPS corresponding to the guess must be left-canonized (+1),
% while the output MPS will be right-canonized (-1).
%
% INPUT	
%	mps_in:				cell array corresponding to input MPS
%	mpo:				cell array corresponding to MPO
%	mps_out:			cell array corresponding to guess of output MPS
%						(WARNING: must be left canonized!)
%	iter_max:			(optional) maximum number of iterations
%	tolerance:			(optional) stopping condition on the relative 
%						improvement of the fidelity. The iterative routine
%						is stopped at the k-th iteration when
%						|fidelity[k] - fidelity[k-1]| < tolerance*fidelity[k]
%	max_storage_size:	(optional) maximum number of bytes available to 
%						store the full application of the MPO to  the MPS.
%						If the actual size is smaller that this limit, the
%						product MPO*MPS will be cached
% OUTPUT
%	mps_out:			resulting MPS after computation in right canonization
%	mps_norm:			overlap between the compressed MPS and the full MPS,
%						for normalized states, the compression error can be 
%						estimated as err = 1 - mps_norm
%	iter:				number of iterations in the optimization

% Handle optional arguments
if nargin < 6
	max_storage_size = inf;
end
if nargin < 5
	tolerance = 1e-6;
end
if nargin < 4
	iter_max = 100;
end

if iter_max <= 0
	iter = 0;
	return
end

N = length(mps_in);

if isempty(mpo)
	mult = @(W,M,k) M;
	% Create dummy MPO
	mpo = cell(1,N);
else
	if predict_product_size(mpo,mps_in) <= max_storage_size
		mps_full = cellfun(@(W,M) apply(W,M),mpo,mps_in,'UniformOutput',false);
		mult = @(W,M,k) mps_full{k};
	else
		mult = @(W,M,k) apply(W,M);
	end
end

% Initialize block storage
blocks = cell(1,N+1);
blocks{1} = 1;
blocks{N+1} = 1;
for site = N:(-1):2
	target = mult(mpo{site},mps_in{site},site);
	blocks{site} = update_block(blocks{site+1},mps_out{site},[],target,-1);
end

mps_norm_prev = 0;
for iter = 1:iter_max
	% Optimization sweep left -> right
	for site = 1:(N-1)
		% Apply MPO
		target = mult(mpo{site},mps_in{site},site);
		% Variational step
		mps_out{site} = optimization_step(target,blocks{site},blocks{site+1});
		% Canonize the new element
		mps_out{site} = canonize_fast(mps_out{site},+1);
		% Update current block
		blocks{site+1} = update_block(blocks{site},mps_out{site},[],target,+1);
	end
	% Do same for last site, except update
	target = mult(mpo{N},mps_in{N},N);
	mps_out{N} = optimization_step(target,blocks{N},blocks{N+1});
	mps_out{N} = canonize_fast(mps_out{N},+1);
	% Optimization sweep right -> left
	for site = N:(-1):2
		% Apply MPO
		target = mult(mpo{site},mps_in{site},site);
		% Optimization step
		mps_out{site} = optimization_step(target,blocks{site},blocks{site+1});
		% Canonize the new element
		mps_out{site} = canonize_fast(mps_out{site},-1);
		% Update current block
		blocks{site} = update_block(blocks{site+1},mps_out{site},[],target,-1);

	end
	% Do same for first site, and compute norm
	target = mult(mpo{1},mps_in{1},1);
	mps_out{1} = optimization_step(target,blocks{1},blocks{2});
	mps_out{1} = canonize_fast(mps_out{1},-1);
	carryover = update_block(blocks{2},mps_out{1},[],target,-1);
	mps_norm = abs(carryover);
	mps_out{1} = sign(carryover)*mps_out{1};
	% Calculate stopping condition
	if abs(mps_norm-mps_norm_prev) < tolerance*mps_norm
		break;
	end
	mps_norm_prev = mps_norm;
end
end

function new_guess = optimization_step(target,block_left,block_right)
	new_guess = contract(block_right,2,2,target,3,2);
	new_guess = contract(block_left,2,2,new_guess,3,2);
end

function bytes = predict_product_size(mpo,mps)
% Predicts the size of the resulting MPS in bytes if you compute the
% application of some MPO to an MPS provided. Notice that this is a 
% slightly optimistic estimate since it doesn't take in account the 
% MATLAB overhead.
% 
% INPUT
%	O:		1d cell array representing an MPO 
%	M:		1d cell array representing an MPS
% OUTPUT
%	nBytes:	number of bytes corresponding to MPO*MPS

complexToBytes = 16;

nComplex = sum(cellfun(@(W,M) size(W,1)*size(W,2)*prod(size(M)),mpo,mps));
bytes = complexToBytes*nComplex;
end
