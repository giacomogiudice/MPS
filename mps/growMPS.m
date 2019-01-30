function [mps,energy,iter] = growMPS(mps,mpo,D_new)
% Handle optional arguments
N = length(mps);

% Initialize block storage
blocks = cell(1,N+1);
blocks{1} = 1;
blocks{N+1} = 1;
for site = N:(-1):2
	blocks{site} = update_block(blocks{site+1},mps{site},mpo{site},mps{site},-1);
end
% Sweep left -> right
for site = 1:(N-1)
	[D1,D2,d] = size(mps{site});
	% If maximum size allowed by OBC, canonize and move on to next site
	D_max = min(d^site,d^(N-site));
	if D2 == D_max
		[mps{site},carryover] = canonize_fast(mps{site},1);
		mps{site+1} = ncon({mps{site+1},carryover},{[1,-2,-3],[-1,1]});
		continue
	end
	% Get nullspaces
	[N_left,mps_left,carryover] = nullspace(mps{site},1);
	N_right = nullspace(mps{site+1},-1);
	% Left and right environment
	G_left = update_block(blocks{site},N_left,mpo{site},mps{site},+1);
	G_right = update_block(blocks{site+2},N_right,mpo{site+1},mps{site+1},-1);
	% Contract environments and compute SVD 
	M = ncon({G_left,G_right},{[-1,1,2],[-2,1,2]});
	[U,~,~] = svd(M);
	D3 = size(mps{site+1},2);
	D_diff = min([D_new,size(U,1)+D2,D_max]) - D2;
	if D_diff < 0
		error('Cannot reduce bond dimension.')
	end
	U = U(:,1:D_diff);
	% Contract carryover with next site
	mps_next = ncon({mps{site+1},carryover},{[1,-2,-3],[-1,1]});
	% Increase bond dimension
	mps{site} = zeros([D1,D2+D_diff,d]);
	mps{site+1} = zeros([D2+D_diff,D3,d]);
	for s = 1:d
		mps{site}(:,:,s) = [mps_left(:,:,s),N_left(:,:,s)*U];
		mps{site+1}(:,:,s) = [mps_next(:,:,s);zeros([D_diff,D3])];
	end
	% Update block
	blocks{site+1} = update_block(blocks{site},mps{site},mpo{site},mps{site},+1);
end
% Left-canonize last site
mps{N} = canonize_fast(mps{N},1);
end
