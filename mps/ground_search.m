function [mps,energy,iter] = ground_search(mps,mpo,iter_max,precision,verbose)
% In the spirit of DMRG, finds an MPS approximation of the state which 
% minimizes of the function <mps|mpo|mps>/<mps|mps>. This problem is
% iteratively approximated by solving a local eigenvalue problem.
% The algorithm should converge monotonically towards a local minima,
% convergence towards a global minimum is not guaranteed. The result 
% should be checked with different randomized initial states and/or larger
% bond dimension. If the algorithm has difficulties converging, the 
% precision might be too demanding or the bond dimension might be too low.
% WARNING: The input MPS must be left-canonized (+1), while the output MPS
% will be right-canonized (-1).
%
% INPUT
%	mps_in:		cell array corresponding to input MPS, usually a random 
%				MPS with sufficient bond dimension is a good start
%				(WARNING: must be left canonized!)
%	mpo:		cell array corresponding to MPO, assumed to be Hermitian           
%	iter_max:	(optional) maximum number of iterations
%	precision:	(optional) stopping condition on the relative improvement
%				of the eigenvalue. The iterative routine is stopped at the
%				k-th iteration when
%				|energy[k] - energy[k-1]| < precision*max(|energy[k]|,1);
%	verbose:	(optional) setting this to true will output on the screen
%				the results at each iteration
% OUTPUT
%	mps:		approximation of the eigenstate in right canonization
%	energy:		approximation of the smallest real eigenvalue		
%	iter:		number of iterations used in the optimization

% Handle optional arguments
if nargin < 5
	verbose = false;
end
if nargin < 4
	precision = eps;
end
if nargin < 3
	iter_max = 100;
end

N = length(mps);

% Initialize block storage
blocks = cell(1,N+1);
blocks{1} = 1;
blocks{N+1} = 1;
for site = N:(-1):2
	blocks{site} = update_block(blocks{site+1},mps{site},mpo{site},mps{site},-1);
end
energy_prev = update_block(blocks{2},mps{1},mpo{1},mps{1},-1);
if verbose
	fprintf('Iter\t      Energy\t Energy Diff\tLap Time [s]\n')
end
for iter = 1:iter_max
	if verbose, tic; end
	% Sweep left -> right
	for site = 1:(N-1)
		% Compute 'one-site' effective Hamiltonian
		fun = fun_onesite(mpo{site},blocks{site},blocks{site+1});
		% Solve the local problem
		mps{site} = optimization_step(mps{site},fun);
		% Canonize the new element
		mps{site} = canonize_fast(mps{site},+1);
		% Compute block update
		blocks{site+1} = update_block(blocks{site},mps{site},mpo{site},mps{site},+1);
	end

	% Sweep right -> left
	for site = N:(-1):2
		% Compute 'one-site' effective Hamiltonian
		fun = fun_onesite(mpo{site},blocks{site},blocks{site+1});
		% Solve the local problem
		[mps{site},energy] = optimization_step(mps{site},fun);
		% Canonize the new element
		[mps{site},carryover] = canonize_fast(mps{site},-1);
		% Compute block update
		blocks{site} = update_block(blocks{site+1},mps{site},mpo{site},mps{site},-1);
	end
	% Print information to screen
	if verbose
		fprintf('%4d\t%12g\t%12g\t%12.1f\n',iter,energy,energy_prev - energy,toc);
	end
	% Stopping condition on the improvement
	if iter > 1 && abs(energy - energy_prev) < precision*max(abs(energy),1)
		break;
	end
	energy_prev = energy;
end
% Right-canonize last site
mps{1} = ncon({mps{1},carryover},{[-1,1,-3],[1,-2]});
mps{1} = canonize_fast(mps{1},-1);
end
