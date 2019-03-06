% This script runs the time-evolution of the Heisenberg model with 
% conventional matrix-vector operations. This is the only routine that does
% not use Matrix Product States (MPS). Because of the exponential increase
% in Hilbert space this method is limited to small (<20) systems.
%% Wipe Everything And Get Parameters
clear
setup

%% Building the Conventional Hamiltonian
H = sparse(0);
for site_1 = 1:N-1
	O_x = 1;
	O_y = 1;
	O_z = 1;
	for site_2 = 1:N
		if site_2 == site_1 || site_2 == site_1 + 1
			O_x = kron(O_x,sigma.x);
			O_y = kron(O_y,sigma.y);
			O_z = kron(O_z,sigma.z);
		else
			O_x = kron(O_x,sigma.id);
			O_y = kron(O_y,sigma.id);
			O_z = kron(O_z,sigma.id);
		end
	end
	H = H + O_x + O_y + O_z ;
end

%% Build the Initial State
state = [1;0];
for site = 2:N
	state = kron(state,[0;1]);
end

%% Build The Magnetization Observable
magnetization = cell(N,1);
for site = 1:N
	magnetization{site} = kron(speye(2^(site-1)),sigma.z);
	magnetization{site} = kron(magnetization{site},speye(2^(N-site)));
end

%% Conventional Evolution
% Array in which to store observables
magn_conventional = zeros(N,time_steps);
% Values for t=0
for site = 1:N
	magn_conventional(site,1) = state'*magnetization{site}*state;
end

% Arnoldi exponential approximation
for step = 2:time_steps
	state = exp_arnoldi(state,@(v) -1i*(H*v),dt);
	state = state/norm(state);
	for site = 1:N
		magn_conventional(site,step) = state'*magnetization{site}*state;
	end
end

%% Save To File
save(filename,'magn_conventional','-append');
