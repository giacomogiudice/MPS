function [mps,energy,iter] = projected_search(mps,mpo,mps_list,iter_max,precision,verbose)
% Performs a minimization of <mps|mpo|mps>/<mps|mps> with the constraint(s)
% <mps|mps_list{n}> = 0 for every element in mps_list. In order to ensure
% that the optimizer actually find a minimum of the energy function, one
% should make sure that this minimum is well below zero, i.e. by adding a
% negative constant to the Hamiltonian.
% WARNING: The input MPS must be left-canonized (+1), as well as the 
% elements in mps_list. The output MPS will be right-canonized (-1).
%
% INPUT
%	mps_in:		cell array corresponding to input MPS, usually a random 
%				MPS with sufficient bond dimension is a good start
%				(WARNING: must be left canonized!)
%	mpo:		cell array corresponding to MPO, assumed to be Hermitian
%	mps_list:	cell array of MPS that correspond to the orthogonal constraints 
%				(WARNING: all elements must be left canonized!)        
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
if nargin < 6
	verbose = false;
end
if nargin < 5
	precision = eps;
end
if nargin < 4
	iter_max = 100;
end

N = length(mps);
M = length(mps_list);
if M == 0
	[mps,energy,iter] = ground_search(mps,mpo,iter_max,precision,verbose);
	return
end
if M > 1 & size(mps_list) == 1
	mps_list = transpose(mps_list);
end
mps_list = cat(1,mps_list{:});

% Initialize block storage
Hblocks = cell(1,N+1);
Hblocks{1} = 1;
Hblocks{N+1} = 1;
for site = N:(-1):2
	Hblocks{site} = update_block(Hblocks{site+1},mps{site},mpo{site},mps{site},-1);
end
energy_prev = update_block(Hblocks{2},mps{1},mpo{1},mps{1},-1);

Cblocks = cell(M,N+1);
for m = 1:M
	Cblocks{m,1} = 1;
	Cblocks{m,N+1} = 1;
	for site = N:(-1):2
		Cblocks{m,site} = update_block(Cblocks{m,site+1},mps{site},{},mps_list{m,site},-1);
	end
end

if verbose
	fprintf('Iter\t      Energy\t Energy Diff\tLap Time [s]\n')
end
for iter = 1:iter_max
	if verbose, tic; end
	% Sweep left -> right
	for site = 1:(N-1)
		% Compute 'one-site' effective Hamiltonian
		ham = fun_onesite(mpo{site},Hblocks{site},Hblocks{site+1});
		proj = proj_onesite(mps_list(:,site),Cblocks(:,site),Cblocks(:,site+1));
		fun = @(M) proj(ham(proj(M)));
		% Solve the local problem
		mps{site} = optimization_step(mps{site},fun);
		% Canonize the new element
		mps{site} = canonize_fast(mps{site},+1);
		% Compute block update
		Hblocks{site+1} = update_block(Hblocks{site},mps{site},mpo{site},mps{site},+1);
		for m = 1:M
			Cblocks{m,site+1} = update_block(Cblocks{m,site},mps{site},{},mps_list{m,site},+1);
		end
	end

	% Sweep right -> left
	for site = N:(-1):2
		% Compute 'one-site' effective Hamiltonian
		ham = fun_onesite(mpo{site},Hblocks{site},Hblocks{site+1});
		proj = proj_onesite(mps_list(:,site),Cblocks(:,site),Cblocks(:,site+1));
		fun = @(M) proj(ham(proj(M)));
		% Solve the local problem
		[mps{site},energy] = optimization_step(mps{site},fun);
		% Canonize the new element
		[mps{site},carryover] = canonize_fast(mps{site},-1);
		% Compute block update
		Hblocks{site} = update_block(Hblocks{site+1},mps{site},mpo{site},mps{site},-1);
		for m = 1:M
			Cblocks{m,site} = update_block(Cblocks{m,site+1},mps{site},{},mps_list{m,site},-1);
		end
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
