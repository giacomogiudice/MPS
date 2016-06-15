function [U_odd,U_even] = evolution_heisenberg(N,J,U,dt)
% Returns odd and even evolution operators for the Heisenberg
% Hamiltonian of parameters N,J,U,dt
d = 2;

%% Pauli Matrices
S_x =[0,1;1,0];
S_y = [0,-1i;1i,0];
S_z = [1,0;0,-1];
%% Pair Evolution Operator
Ham_pair = J*(kron(S_x,S_x) + kron(S_y,S_y)) + U*(kron(S_z,S_z));
U_pair = expm(-1i*dt*Ham_pair);

%% U -> P reshaping
P = reshape(U_pair,[2,2,2,2]);
P = permute(P,[4,2,3,1]); 
P = reshape(P,[d^2,d^2]); 

%% SVD of P
[U,S,V] = svd(P,'econ');
k = size(S,2);
U = U*sqrt(S); 
U_bar = sqrt(S)*V'; 

%% U, U_bar and eye(2) reshaping
U = reshape(U,[d,d,1,k]);
U = permute(U,[3,4,1,2]);
U_bar = reshape(U_bar,[k,1,d,d]);
I_site = reshape(eye(2),[1,1,d,d]);

%% Putting it all into U_even and U_odd
U_odd = cell(1,N);
for i = 1:2:N-1;
U_odd{i}= U;
U_odd{i+1}= U_bar;
end

U_even = cell(1,N);
U_even{1}= I_site;
for i = 2:2:N-1;
U_even{i} = U;
U_even{i+1} = U_bar;
end

% Fix chain length problems
if not(mod(N,2))
    U_even{N} = I_site;
else
    U_odd{N} = I_site;
end

%% Compare with Trotter Function
[U_o, U_e] = trotter(U_pair,N,d);
err = 0;
for i = 1:N
    err = err + norm(U_o{i}(:) - U_odd{i}(:))+ norm(U_e{i}(:) - U_even{i}(:));
end
assert(err < 1e-12);
end