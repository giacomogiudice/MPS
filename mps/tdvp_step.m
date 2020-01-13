function [mps,blocks,energy] = tdvp_step(mps,mpo,dt,blocks,settings)

N = length(mps);
dt_half = dt/2;
integrator = settings.integrator.handle;
integrator_options = settings.integrator.options;

for site = 1:(N-1)
    % Compute 'one-site' effective Hamiltonian
    fun = fun_onesite(mpo{site},blocks{site},blocks{site+1});
    % Evolve the current element forward
    mps{site} = integrator(dt_half,fun,mps{site},integrator_options);
    % Canonize the new element
    [mps{site},carryover] = canonize_fast(mps{site},+1);
    % Compute block update
    new_block = update_block(blocks{site},mps{site},mpo{site},mps{site},+1);
    % Compute 'zero-site' effective Hamiltonian
    fun = fun_zerosite(new_block,blocks{site+1});
    % Evolve the carryover backward
    carryover = integrator(-dt_half,fun,carryover,integrator_options);
    % Finally contract the carryover with the next site and update block
    mps{site+1} = ncon({carryover,mps{site+1}},{[-1,1],[1,-2,-3]});
    blocks{site+1} = new_block;
end
% Do only the forward step for the last site
fun = fun_onesite(mpo{N},blocks{N},blocks{N+1});
mps{N} = integrator(dt_half,fun,mps{N},integrator_options);
[mps{N},carryover] = canonize_fast(mps{N},+1);

% Sweep right -> left
for site = N:(-1):2
    % Compute 'one-site' effective Hamiltonian
    fun = fun_onesite(mpo{site},blocks{site},blocks{site+1});
    % Evolve the current element forward
    mps{site} = integrator(dt_half,fun,mps{site},integrator_options);
    % Canonize the new element
    [mps{site},carryover] = canonize_fast(mps{site},-1);
    % Compute block update
    new_block = update_block(blocks{site+1},mps{site},mpo{site},mps{site},-1);
    % Compute 'zero-site' effective Hamiltonian
    fun = fun_zerosite(blocks{site},new_block);
    % Evolve the carryover backward
    carryover = integrator(-dt_half,fun,carryover,integrator_options);
    % Finally contract the carryover with the next site and update block
    mps{site-1} = ncon({mps{site-1},carryover},{[-1,1,-3],[1,-2]});
    blocks{site} = new_block;
end
% Do only the forward step for the first site
fun = fun_onesite(mpo{1},blocks{1},blocks{2});
mps{1} = integrator(dt_half,fun,mps{1},integrator_options);
[mps{1},carryover] = canonize_fast(mps{1},-1);
energy = update_block(blocks{2},mps{1},mpo{1},mps{1},-1);
end

