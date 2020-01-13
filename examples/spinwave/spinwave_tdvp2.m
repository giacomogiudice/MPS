% This script runs the time-evolution of the Heisenberg model using
% the Time-Dependent Variational Principle (TDVP).
%% Wipe Everything And Get Parameters
clear
setup

%% Define Hamiltonian Operator
D_O = 5;
W = zeros(D_O,D_O,d,d);
W(1,1,:,:) = sigma.id;
W(2,1,:,:) = sigma.x;
W(3,1,:,:) = sigma.y;
W(4,1,:,:) = sigma.z;

W(5,2,:,:) = sigma.x;
W(5,3,:,:) = sigma.y;
W(5,4,:,:) = sigma.z;
W(5,5,:,:) = sigma.id;

H = cell(1,N);
H{1} = W(D_O,:,:,:);
for i = 2:N-1
    H{i} = W;
end
H{N} = W(:,1,:,:);

%% Build Initial State
state = cell(1,N);
for site = 1:N
    state{site} = reshape([0,1],[1,1,d]);
end
state{1} = reshape([1,0],[1,1,d]);
state{2} = reshape([1,0],[1,1,d]);
% Pad state
state = sweep(state,{},-1);
% state = padMPS(state,D_static,-1);

%% Time Evolution
magn_tdvp2 = zeros(N,time_steps);
% Values for t=0
magn_tdvp2(:,1) = real(expectationvalue(state,sigma.z));
% Trotter-Suzuki order 2
for step = 2:time_steps
    state = tdvp2_step(state,H,-1i*dt,D_static,eps);
    state
    magn_tdvp2(:,step) = real(expectationvalue(state,sigma.z));
end

%% Save To File
save(filename,'magn_tdvp2','-append');
