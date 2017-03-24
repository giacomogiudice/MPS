function mps = tdvp_step(mps,mpo,dt)
% blah

N = length(mps);
dt_half = dt/2;

% Initialize block storage
blocks = cell(1,N+1);
blocks{1} = 1;
blocks{N+1} = 1;
for site = N:(-1):2
	blocks{site} = update_block(blocks{site+1},mps{site},mpo{site},mps{site},-1);
end

% Sweep left -> right
for site = 1:(N-1)
	% Compute 'one-site' effective Hamiltonian
	fun = ham_onesite(mpo{site},blocks{site},blocks{site+1});
	% Evolve the current element forward
	mps{site} = RK4_step(mps{site},fun,dt_half);
    % Canonize the new element
    [mps{site},carryover] = canonize_fast(mps{site},+1);
	% Compute block update
	new_block = update_block(blocks{site},mps{site},mpo{site},mps{site},+1);
	% Compute 'zero-site' effective Hamiltonian
	fun = ham_zerosite(new_block,blocks{site+1});
	% Evolve the carryover backward
	carryover = RK4_step(carryover,fun,dt_half);
	% Finally contract the carryover with the next site and update block
	mps{site+1} = contract(carryover,2,2,mps{site+1},3,1);
	blocks{site+1} = new_block;
end
% Do only the forward step for the last site
fun = ham_onesite(mpo{N},blocks{N},blocks{N+1});
mps{N} = RK4_step(mps{N},fun,dt_half);
mps{N} = canonize_fast(mps{N},+1);

% Sweep right -> left
for site = N:(-1):2
	% Compute 'one-site' effective Hamiltonian
	fun = ham_onesite(mpo{site},blocks{site},blocks{site+1});
	% Evolve the current element forward
	mps{site} = RK4_step(mps{site},fun,dt_half);
    % Canonize the new element
    [mps{site},carryover] = canonize_fast(mps{site},-1);
	% Compute block update
	new_block = update_block(blocks{site+1},mps{site},mpo{site},mps{site},-1);
	% Compute 'zero-site' effective Hamiltonian
	fun = ham_zerosite(blocks{site},new_block);
	% Evolve the carryover backward
	carryover = RK4_step(carryover,fun,dt_half);
	% Finally contract the carryover with the next site and update block
	mps{site-1} = permute(contract(mps{site-1},3,2,carryover,2,1),[1 3 2]);
	blocks{site} = new_block;
end
% Do only the forward step for the first site
fun = ham_onesite(mpo{1},blocks{1},blocks{2});
mps{1} = RK4_step(mps{1},fun,dt_half);
mps{1} = canonize_fast(mps{1},-1);

end

function handle = ham_onesite(op,left,right)
	function W = new_element(M)
		W = contract(right,3,2,M,3,2);
		W = contract(W,4,[2,4],op,4,[2,4]);
		W = contract(left,3,[2,3],W,4,[2,3]);
		W = -1i*W;
	end
	handle = @new_element;
end

function handle = ham_zerosite(left,right)
	function W = new_element(C)
		W = contract(right,3,2,C,2,2);
		W = contract(left,3,[2,3],W,3,[3,2]);
		W = 1i*W;
	end
	handle = @new_element;
end

function v_out = RK4_step(v_in,fun,dt)
	k1 = fun(v_in);
	k2 = fun(v_in + (dt/2)*k1);
	k3 = fun(v_in + (dt/2)*k2);
	k4 = fun(v_in + dt*k3);
	v_out = v_in + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
end
