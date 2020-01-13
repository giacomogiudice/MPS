function [mps,blocks,energy] = tdvp2_step(mps,mpo,dt,blocks,settings)
N = length(mps);
dt_half = dt/2;
integrator = settings.integrator.handle;
integrator_options = settings.integrator.options;
epsilon = settings.truncation;
D_max = settings.maxbond;

% Sweep left -> right
for site = 1:(N-2)
    % Get current 'two-site' block
    theta = merge_twosite(mps{site:site+1});
    % Compute 'two-site' effective Hamiltonian
    fun = fun_twosite(mpo(site:site+1),blocks{site},blocks{site+2});
    % Evolve the current element forward
    theta = integrator(dt_half,fun,theta,integrator_options);
    % Split 'two-site' block and get next tensor
    [mps{site},schmidt,carryover] = split_twosite(theta,D_max,epsilon);
    carryover = ncon({schmidt,carryover},{[-1,1],[1,-2,-3]});
    % Compute block update
    new_block = update_block(blocks{site},mps{site},mpo{site},mps{site},+1);
    % Compute 'one-site' effective Hamiltonian
    fun = fun_onesite(mpo{site+1},new_block,blocks{site+2});
    % Evolve the carryover backward and update block
    mps{site+1} = integrator(-dt_half,fun,carryover,integrator_options);
    blocks{site+1} = new_block;
end
% Do only the forward step for the last two sites
theta = merge_twosite(mps{N-1:N});
 % Compute 'two-site' effective Hamiltonian
fun = fun_twosite(mpo(N-1:N),blocks{N-1},blocks{N+1});
% Evolve the current element forward
theta = integrator(dt_half,fun,theta,integrator_options);
% Split 'two-site' block and get both tensors
[mps{N-1},schmidt,carryover] = split_twosite(theta,D_max,epsilon);
mps{N} = ncon({schmidt,carryover},{[-1,1],[1,-2,-3]});

% Sweep right -> left
for site = (N-1):(-1):2
    % Get current 'two-site' block
    theta = merge_twosite(mps{site:site+1});
    % Compute 'two-site' effective Hamiltonian
    fun = fun_twosite(mpo(site:site+1),blocks{site},blocks{site+2});
    % Evolve the current element forward
    theta = integrator(dt_half,fun,theta,integrator_options);
    % Split 2-site block and get next tensor
    [carryover,schmidt,mps{site+1}] = split_twosite(theta,D_max,epsilon);
    carryover = ncon({carryover,schmidt},{[-1,1,-3],[1,-2]});
    % Compute block update
    new_block = update_block(blocks{site+2},mps{site+1},mpo{site+1},mps{site+1},-1);
    % Compute 'one-site' effective Hamiltonian
    fun = fun_onesite(mpo{site},blocks{site},new_block);
    % Evolve the carryover backward and update block
    mps{site} = integrator(-dt_half,fun,carryover,integrator_options);
    blocks{site+1} = new_block;
end
% Do only the forward step for the first two site
theta = merge_twosite(mps{1:2});
 % Compute 'two-site' effective Hamiltonian
fun = fun_twosite(mpo(1:2),blocks{1},blocks{3});
% Evolve the current element forward
theta = integrator(dt_half,fun,theta,integrator_options);
% Split 'two-site' block and get both tensors
[carryover,schmidt,mps{2}] = split_twosite(theta,D_max,epsilon);
mps{1} = canonize_fast(ncon({carryover,schmidt},{[-1,1,-3],[1,-2]}),-1);
% Compute block update
blocks{2} = update_block(blocks{3},mps{2},mpo{2},mps{2},-1);
energy = update_block(blocks{2},mps{1},mpo{1},mps{1},-1);
end


function theta = merge_twosite(M1,M2)
theta = ncon({M1,M2},{[-1,1,-2],[1,-4,-3]});
end

function [A_left,S,A_right] = split_twosite(theta,D_max,epsilon)
[D_1,d_1,d_2,D_2] = size(theta);
C = reshape(theta,[D_1*d_1,d_2*D_2]);
[U,S,V] = svd(C,'econ');
[U,S,V,D_cut] = trim(U,S,V,D_max,epsilon);
A_left = permute(reshape(U,[D_1,d_1,D_cut]),[1,3,2]);
A_right = permute(reshape(V',[D_cut,d_2,D_2]),[1,3,2]);
end


function [U,S,V,D_cut] = trim(U,S,V,D_max,epsilon)
% Computes appropriate trimming of U,S,V matrices after SVD
% The decimation is computed by taking the sum of the square of the singular
% values and then keeping only the ones under (1 - epsilon)
%
% INPUT
%   U,S,V:      matrices from SVD decompositon
%   D_max:      maximum allowed bond length
%   epsilon:    tolerated error in decimation
% OUTPUT
%   U,S,V:      matrices from SVD decompositon after truncation
%   D_cut:      computed bond length

% Compute number of singular values under (1 - epsilon)
sum_s = cumsum(diag(S).^2);
D_cut = sum(sum_s < (1 - epsilon)*sum_s(end)) + 1;
if D_cut > D_max
    D_cut = D_max;
end
% Do trimming
U = U(:,1:D_cut);
S = S(1:D_cut,1:D_cut);
V = V(:,1:D_cut);
end
