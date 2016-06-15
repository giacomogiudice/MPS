% This stript test 'trotter' against the 'evolution_heisenberg'

close all
clear 
%% Initial Parameters

N = 8;
J = 1;
U = 1;
d = 2;
dt = 0.1; % Time increment

[U_odd,U_even] = evolution_heisenberg(N,J,U,dt);
%% Stuff
S_X =[0,1;1,0];
S_Y = [0,-1i;1i,0];
S_Z = [1,0;0,-1];
Ham_pair = J*(kron(S_X,S_X) + kron(S_Y,S_Y)) + U*(kron(S_Z,S_Z));
U_pair = expm(-1i*dt*Ham_pair);
%% Kron
 H_odd_kron = 0;

 for i = 1:2:N-1
    H_odd_kron = H_odd_kron + kron(kron(eye(d^(i-1)),Ham_pair),eye(d^(N-i-1)));
 end
 U_odd_kron = expm(-1i*dt*H_odd_kron);

 H_even_kron = 0;
for i = 2:2:N-1
    H_even_kron = H_even_kron + kron(kron(eye(d^(i-1)),Ham_pair),eye(d^(N-i-1)));
end
U_even_kron = expm(-1i*dt*H_even_kron);
%% MPS U
fprintf('Testing U_odd...');
assert(norm(expand_mpo(U_odd) - U_odd_kron) < 1e-12);
fprintf('\tdone\n');
fprintf('Testing U_even...');
assert(norm(expand_mpo(U_even) - U_even_kron) < 1e-12);
fprintf('\tdone\n');
[U_odd_t,U_even_t] = trotter(U_pair,N,d);
fprintf('Comparing with "trotter.m" ...');
assert(norm(expand_mpo(U_odd) - expand_mpo(U_odd_t)) + norm(expand_mpo(U_even) - expand_mpo(U_even_t)) < 1e-12);
fprintf('\tdone\n');