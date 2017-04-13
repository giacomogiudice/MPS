% This script runs the time-evolution of the Heisenberg model with 
% Trotter decomposition and then iterative compression using a random
% initial guess. One can experiment with better guesses, such as using
% the previous time step (correctly padded to allow the bond dimension to 
% grow) or the SVD compression of MPO*MPS (which is an inexact, but often
% conservative truncation, which provides an excellent initial guess).
%% Wipe Everything And Get Parameters
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
magn_iter = zeros(N,time_steps);
compression_iter = zeros(1,time_steps);
compression_err = zeros(1,time_steps);

% Values for t=0
state = sweep(state,{},-1);
magn_iter(:,1) = real(expectationvalue(magnetization,state));

% Trotter-Suzuki order 2
for step = 2:time_steps
	[state,iter,err] = sweep_iter(state,U,randomMPS(N,D_static,d,1),iter_max,tolerance);
	magn_iter(:,step) = real(expectationvalue(magnetization,state));
	compression_iter(step) = iter;
	compression_err(step) = err;
end

%% Save To File
save(filename,'magn_iter','compression_iter','compression_err','-append');
