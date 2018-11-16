function g = error_variance(mps,mpo)
% Assumes right canonization (-1)
N = length(mps);
d = size(mps{1},3);

% Initialize blocks storage
blocks = cell(1,N+1);
blocks{1} = 1;
blocks{N+1} = 1;
for site = N:(-1):2
	blocks{site} = update_block(blocks{site+1},mps{site},mpo{site},mps{site},-1);
end

g = 0;
for site = 1:(N-1)
	N_left = nullspace(mps{site},1);
	G_left = update_block(blocks{site},N_left,mpo{site},mps{site},+1);
	% First projector
	P1 = ncon({G_left,blocks{site+1}},{[-1,1,2],[-2,1,2]});
	[mps{site},carryover] = canonize_fast(mps{site},+1);
	blocks{site+1} = update_block(blocks{site},mps{site},mpo{site},mps{site},+1);
	N_right = nullspace(mps{site+1},-1);
	G_right = update_block(blocks{site+2},N_right,mpo{site+1},mps{site+1},-1);
	% Second projector
	P2 = ncon({G_left,G_right},{[-1,1,2],[-2,1,2]});
	g = g +  trace(P1*P1') + trace(P2*P2');
	mps{site+1} = ncon({mps{site+1},carryover},{[1,-2,-3],[-1,1]});
end
% First projector only on last site
N_left = nullspace(mps{N},1);
G_left = update_block(blocks{N},N_left,mpo{N},mps{N},+1);
P1 = ncon({G_left,blocks{N+1}},{[-1,1,2],[-2,1,2]});
g = g + trace(P1*P1');
end