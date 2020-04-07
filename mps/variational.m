function [mps,output,mps_left] = variational(mpo,D_list,settings)
% Performs single-site DMRG, by variationally finding an MPS approximation
% of the state which minimizes (or optionally maximizes) of the function
% <mps|mpo|mps>/<mps|mps>. The behavior can be tuned by passing a settings
% struct, see 'variational_settings.m'. In order too change a default
% setting, simply create a structure with the settings you wish to change,
% and pass it as an argument. The algorithm should converge towards a local
% extremum, bu convergence towards a global minimum is not guaranteed.
% If the algorithm has difficulties converging, the precision might be too
% demanding or the bond dimension might be too low.
%
% INPUT
%	mpo:		cell array corresponding to MPO, assumed to be Hermitian
%	D_list:		array corresponding to the bond dimensions to use for the
%				variational MPS. The routine will automatically increase
%				the bond dimension of the guess so that the final MPS will
%				have a bond dimension corresponding to the last entry of
%				this array. The bond dimension is increased at
%				logarithmically-spaced steps between 1 and settings.tol.
%				(WARNING: The array should be in ascending order).
%	settings:	structure corresponding to the options for the routine
% OUTPUT
%	mps:		approximation of the eigenstate in right canonization (-1)
%	output:		approximation of the smallest real eigenvalue
%	mps_left:	approximation of the eigenstate in left canonization (+1)

% Extract settings
if nargin == 2
	settings = variational_settings();
else
	settings = variational_settings(settings);
end
assert(iscell(mpo),'MPO must be a cell array.');
N = length(mpo);
assert(N > 0,'MPO size must be a positive integer.');
d = size(mpo{1},3);

% Do checks on D_list
assert(all(diff(D_list) > 0),'D_list must be a vector of positive integers in ascending order.');
D_now = D_list(1);
bond_ind = 1;
growbond = false;
growtol = logspace(0,log10(settings.tol),length(D_list)+1);
growtol(1) = [];
growtol(end) = 0;

% Check if initial mps is provided
if isfield(settings.initial,'mps')
	mps = settings.initial.mps;
	assert(length(mps) == N,'Size mismatch between input MPO and initial MPS.');
	D_now = max(cellfun(@(M) max(size(M)),mps));
	assert(D_now <= D_list(1),'Bond dimension of initial MPS must be smaller or equal to first element of D_list.');
	if D_now < D_list(1)
		growbond = true;
	end
else
	mps = randomMPS(N,D_list(1),d,-1,settings.isreal);
end


% Initialize block storage
blocks = cell(1,N+1);
blocks{1} = 1;
blocks{N+1} = 1;
for site = N:(-1):2
	blocks{site} = update_block(blocks{site+1},mps{site},mpo{site},mps{site},-1);
end
energy_prev = update_block(blocks{2},mps{1},mpo{1},mps{1},-1);

% Check if orthogonal blocks exist
N_orth = 0;
if ~isempty(settings.orthogonalize)
	mps_list = settings.orthogonalize;
	N_orth = size(mps_list,1);
	orth_blocks = cell(N_orth,N+1);
	for m = 1:N_orth
		orth_blocks{m,1} = 1;
		orth_blocks{m,N+1} = 1;
		for site = N:(-1):2
			orth_blocks{m,site} = update_block(orth_blocks{m,site+1},mps{site},{},mps_list{m,site},-1);
		end
	end
end
if settings.verbose
	fprintf('Iter\t      Energy\t Energy Diff\tLap Time [s]\tBond Dim\n')
end
err = [];
for iter = 1:settings.maxit
	if settings.verbose, tic; end
	% Sweep left -> right
	for site = 1:(N-1)
		% Compute 'one-site' effective Hamiltonian
		ham = fun_onesite(mpo{site},blocks{site},blocks{site+1});
		if N_orth
			proj = proj_onesite(mps_list(:,site),orth_blocks(:,site),orth_blocks(:,site+1));
			fun = @(M) proj(ham(proj(M)));
		else
			fun = ham;
		end
		% Solve the local problem
		mps{site} = optimization_step(mps{site},fun,settings);
		% Canonize the new element (possibly increase bond dimension)
		[D_past,D_curr,~] = size(mps{site});
		D_max = min(d^site,d^(N-site));
		if growbond && D_curr < D_max
			% Get nullspaces
			[N_left,mps_left] = nullspace(mps{site},+1);
			N_right = nullspace(mps{site+1},-1);
			% Left and right environment
			G_left = update_block(blocks{site},N_left,mpo{site},mps{site},+1);
			G_right = update_block(blocks{site+2},N_right,mpo{site+1},mps{site+1},-1);
			% Contract environments and compute SVD
			M = ncon({G_left,G_right},{[-1,1,2],[-2,1,2]});
			[U,~,~] = svd(M);
			D_next = size(mps{site+1},2);
			D_diff = min([D_list(bond_ind),size(U,1)+D_curr,D_max]) - D_curr;
			U = U(:,1:D_diff);
			mps_next = mps{site+1};
			mps{site} = zeros([D_past,D_curr+D_diff,d]);
			mps{site+1} = zeros([D_curr+D_diff,D_next,d]);
			for s = 1:d
				mps{site}(:,:,s) = [mps_left(:,:,s),N_left(:,:,s)*U];
				mps{site+1}(:,:,s) = [mps_next(:,:,s);zeros([D_diff,D_next])];
			end
		else
			mps{site} = canonize_fast(mps{site},+1);
		end
		% Compute block update
		blocks{site+1} = update_block(blocks{site},mps{site},mpo{site},mps{site},+1);
		for m = 1:N_orth
			orth_blocks{m,site+1} = update_block(orth_blocks{m,site},mps{site},{},mps_list{m,site},+1);
		end
	end
	% Sweep right -> left
	for site = N:(-1):2
		% Compute 'one-site' effective Hamiltonian
		ham = fun_onesite(mpo{site},blocks{site},blocks{site+1});
		if N_orth
			proj = proj_onesite(mps_list(:,site),orth_blocks(:,site),orth_blocks(:,site+1));
			fun = @(M) proj(ham(proj(M)));
		else
			fun = ham;
		end
		% Solve the local problem
		[mps{site},energy] = optimization_step(mps{site},fun,settings);
		energy = real(energy);
		% Canonize the new element
		[mps{site},carryover] = canonize_fast(mps{site},-1);
		% Compute block update
		blocks{site} = update_block(blocks{site+1},mps{site},mpo{site},mps{site},-1);
		for m = 1:N_orth
			orth_blocks{m,site} = update_block(orth_blocks{m,site+1},mps{site},{},mps_list{m,site},-1);
		end
	end
	err = energy_prev - energy;
	% Print information to screen
	if settings.verbose
		fprintf('%4d\t%12g\t%12g\t%12.1f\t%8d\n',iter,energy,err,toc,D_list(bond_ind));
	end
	growbond = false;
	% Stopping condition on the improvement
	if bond_ind == length(D_list) && abs(err) < settings.tol*max(abs(energy),1)
		break;
	% Check if the bond dimension can be increased
	elseif err < growtol(bond_ind)*max(abs(energy),1) && bond_ind < length(D_list)
		growbond = true;
		bond_ind = bond_ind + 1;
	end
	energy_prev = energy;
end
% Right-canonize last site
mps{1} = ncon({mps{1},carryover},{[-1,1,-3],[1,-2]});
mps{1} = canonize_fast(mps{1},-1);
% Compute variance and left-canonical form
[variance,mps_left] = error_variance(mps,mpo,blocks);

output.iter = iter;
output.err = err;
output.energy = energy;
output.variance = variance;
end
