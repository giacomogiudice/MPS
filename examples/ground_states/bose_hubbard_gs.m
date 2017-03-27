% This script creates the MPO Hamiltonian corresponding to the Bose-Hubbard
% model. The ground state of this system is them found through variational
% optimization. Once this is complete, properties of this ground state are
% plotted, such as the occupation number and the g2 correlation function.
% Try changing J from 0.1 to 0.5 to observe the transition from Mott 
% insulator to superfluid.
%% Wipe Everything And Get Parameters
clear

%% Parameters
N = 20;		% Number of sites
d = 4;		% Truncation of the local Hilbert space (this is enough without changing mu)
D = 20;		% Bond dimension
mu = 0.5;	% Chemical potential, in units of U
J = 0.1;	% Hopping strength, in units of J (transition around 0.2 for mu=0.5)
precision = 1e-8;
max_iter = 15;

% Local operators
id = eye(d);
a = diag(sqrt(1:(d-1)),1);

%% Define Hamiltonian Operator
D_O = 4;
W = zeros(D_O,D_O,d,d);
W(1,1,:,:) = id;
W(2,1,:,:) = a;
W(3,1,:,:) = a';
W(4,1,:,:) = -mu*a'*a + a'*(a'*a)*a;
W(4,2,:,:) = -J*a';
W(4,3,:,:) = -J*a;
W(4,4,:,:) = id;

H = cell(1,N);
H{1} = W(D_O,:,:,:);
for i = 2:N-1
    H{i} = W;
end
H{N} = W(:,1,:,:);

%% Create Observables
identity = cell(1,N);
occupation = cell(N,N);
correlation = cell(N,N);
for site = 1:N
    identity{site} = reshape(id,[1 1 d d]);
end
mid_site = ceil(N/2);
for site = 1:N
	occupation(site,:) = identity;
	occupation{site,site} = reshape(a'*a,[1 1 d d]);
	if site == mid_site
    	correlation{site,mid_site} = reshape(a'*(a'*a)*a,[1 1 d d]);
    else
		correlation{site,mid_site} = reshape(a'*a,[1 1 d d]);
		correlation{site,site} = reshape(a'*a,[1 1 d d]);
	end
end

%% Do Optimization
state = sweep(randomMPS(N,D,d,1),{},-1);
[state,E,iter] = ground_search(state,H,max_iter,precision,true);

%% Plot Observables
n_particles = real(expectationvalue(occupation,state));
g2 = real(expectationvalue(correlation,state))./(n_particles(mid_site).*n_particles);

figure
plot(1:N,n_particles,'s--');
xlabel('$k$')
ylabel('$\langle a^\dagger_k a_k \rangle$')
xlim([0 d]);
figure
hold on
plot(1:N,g2);
xlabel('$k$')
ylabel('$g^2_{{\rm mid},k}$')
