% This script runs the time-evolution of the Heisenberg model with 
% Trotter decomposition and then two sweeps over the state MPS: first
% canonization and  then compression with SVD. In this way the bond 
% dimension is controlled dynamically.
clear
setup

%% Define Evolution Operator
% Two-site Hamiltonian
Ham = kron(sigma.x,sigma.x) + kron(sigma.y,sigma.y) + kron(sigma.z,sigma.z);

U_odd = cell(1,N);
U_even = cell(1,N);
% Perform Trotter decomposition
[V,W] = trotter(expm(-1i*Ham*dt/2),d);
% Build Operators for second-order Trotter-Suzuki
for site = 1:2:N-1
	U_odd{site} = V;
	U_odd{site+1} = W;
end
[V,W] = trotter(expm(-1i*Ham*dt),d);
for site = 2:2:N-1
	U_even{site} = V;
	U_even{site+1} = W;
end
% Fix boundaries
U_even{1} = reshape(sigma.id,[1,1,d,d]);
if mod(N,2)
	U_odd{N} = reshape(sigma.id,[1,1,d,d]);
else
	U_even{N} = reshape(sigma.id,[1,1,d,d]);
end
% Do  compression to merge into a single operator
U = compressMPO(U_odd,U_even,U_odd); 

%% Build Initial State
state = cell(1,N);
for site = 1:N
	state{site} = reshape([0,1],[1,1,d]);
end
state{1} = reshape([1,0],[1,1,d]);

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

maxbond = @(s) max(cellfun(@(x) max(size(x,1),size(x,2)),s));

%% Time Evolution
magn_svd = zeros(N,time_steps);
bond_svd = ones(1,time_steps);
% Values for t=0
magn_svd(:,1) = real(expectationvalue(magnetization,state));
% Trotter-Suzuki order 2
for step = 2:time_steps
	state = sweep(state,U,1);
	state = sweep(state,{},-1,D_max,epsilon);
	magn_svd(:,step) = real(expectationvalue(magnetization,state));
	bond_svd(step) = maxbond(state);
end

%% Save To File
save(filename,'magn_svd','bond_svd','-append');
