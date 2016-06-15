% Another script to test different the time-evolution and block decimation 
% using MPS. It then compares it with the traditional matrix-vector 
% operations. Expected output is a nice little animation comparing the 
%occupation number that evolves in time. 
% Feel free to change all the simulation parameters.

close all
clear
%% Parameters
N = 9;
J = 1;
U = 0;
T = 10;
D_max = 10;
epsilon = 1e-9;
dt = 5e-2;
time = 0:dt:T;

%% Pauli Matrices
S_x = [0,1;1,0];
S_y = [0,-1i;1i,0];
S_z = [1,0;0,-1];
S_plus = [0,1;0,0];
S_minus = S_plus';
I = eye(2);

%% Define MPO Operators
state = cell(1,N);
for i = 1:N
    state{i} = reshape([0,1],[1,1,2]);
end
state{1} = reshape([1,0],[1,1,2]);

[U_odd,U_even] = evolution_heisenberg(N,J,U,dt/2);

I_MPO = cell(1,N);
for i = 1:N
    I_MPO{i} = reshape(I,[1,1,2,2]);
end
S_z_MPO = I_MPO;

%% Do Evolution
profiles = zeros(N,length(time));

j = 1;
for step = 1:length(time)
    for i = 1:N
        S_z_MPO = I_MPO;
        S_z_MPO{i} = reshape(S_z,[1,1,2,2]);
        profiles(i,step) = real(matrixelement(state,S_z_MPO,state));
        S_z_MPO{i} = reshape(eye(2),[1,1,2,2]);
    end
    state = sweep(state,U_odd,+1,D_max,epsilon);
    state = sweep(state,U_even,-1,D_max,epsilon);
    state = sweep(state,U_even,+1,D_max,epsilon);
    state = sweep(state,U_odd,-1,D_max,epsilon);
    j = j+1;
end

%% Building the Conventional Hamiltonian
H = sparse(0);
for i = 1:N-1
    O_x = 1;
    O_y = 1;
    O_z = 1;
    for j = 1:N
        if j == i || j == i+1
            O_x = kron(O_x,S_x);
            O_y = kron(O_y,S_y);
            O_z = kron(O_z,S_x);
        else
            O_x = kron(O_x,I);
            O_y = kron(O_y,I);
            O_z = kron(O_z,I);
        end
    end
    H = H + J*O_x + J*O_y + U*O_z ;
end
S_z_indexed = cell(N,1);
S_number_indexed = cell(N,1);

for i = 1:N
    S_z_i = kron(speye(2^(i-1)),S_z);
    S_z_i = kron(S_z_i,speye(2^(N-i)));
    S_number_i = kron(speye(2^(i-1)),S_plus*S_minus);
    S_number_i = kron(S_number_i,speye(2^(N-i)));
    S_z_indexed{i} = S_z_i;
    S_number_indexed{i} = S_number_i;
end

%% Initial State
state = [1;0];
for k = 2:N
    state = kron(state,[0;1]);
end
%% Conventional Evolution
profiles_conv = zeros(N,length(time));
for step = 1:length(time)
    for i = 1:N
        profiles_conv(i,step) = state'*((S_z_indexed{i})*state);
    end
    %RK-4
    k1 = -1i*(H*state);
    k2 = -1i*H*(state + (dt/2)*k1);
    k3 = -1i*H*(state + (dt/2)*k2);
    k4 = -1i*H*(state + dt*k3);
    state = state + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    state = state/norm(state);
end

%% Animation
figure
sites = 1:N;
for step = 1:length(profiles)
    plot(sites,profiles(:,step),sites,profiles_conv(:,step))
    legend('MPS','conventional');
    axis([0.5 N+0.5 -1.2 1.2])
    pause(1/30)
end
