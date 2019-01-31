% This script creates the MPO Hamiltonian corresponding to the Bose-Hubbard
% model. The ground state of this system is them found through variational
% optimization. Once this is complete, properties of this ground state are
% plotted, such as the occupation number and the correlation correlation function.
% Try changing J from 0.1 to 0.5 to observe the transition from Mott 
% insulator to superfluid.
%% Wipe Everything And Get Parameters
clear

%% Parameters
N = 30;		% Number of sites
d = 4;		% Truncation of the local Hilbert space (this is enough without changing mu)
D = 20;		% Bond dimension
mu = 0.5;	% Chemical potential, in units of U
J = 0.2;	% Hopping strength, in units of J (transition around 0.2 for mu=0.5)

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

%% Do Optimization
settings.tolerance = 1e-8;
settings.maxit = 15;
[state,output] = variational(H,D,settings);

%% Plot Observables
mid_site = ceil(N/2);
n_particles = real(expectationvalue(state,{a'*a}));
correlation = abs(expectationvalue(state,{a',a},[mid_site,N]));

figure(1)
plot(1:N,n_particles,'s--');
xlabel('$k$')
ylabel('$\langle a^\dagger_k a_k \rangle$')
hold on

figure(2)
plot(mid_site:(N-1),correlation);
xlabel('$k$')
ylabel('$|\langle a^\dagger_{\rm mid} a_k|$')
set(gca,'yscale','log')
hold on
