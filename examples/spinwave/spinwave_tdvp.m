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
% Pad state
state = padMPS(state,D_static,-1);

%% Create Observables
identity = cell(1,N);
magnetization = cell(N,N);
for site = 1:N
	identity{site} = reshape(sigma.id,[1 1 d d]);
end
for site = 1:N
	magnetization(site,:) = identity;
	magnetization{site,site} = reshape(sigma.z,[1,1,d,d]);
end

%% Time Evolution
magn_tdvp = zeros(N,time_steps);
% Values for t=0
magn_tdvp(:,1) = real(expectationvalue(magnetization,state));
% Trotter-Suzuki order 2
for step = 2:time_steps
	state = padMPS(state,D_static,-1);
	state = tdvp_step(state,H,dt);
	magn_tdvp(:,step) = real(expectationvalue(magnetization,state));
end

%% Save To File
save(filename,'magn_tdvp','-append');
