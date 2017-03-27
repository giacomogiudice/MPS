% This script creates the MPO Hamiltonian corresponding to the Heisenberg
% model, with an external field of strength g along x. The ground state of
% this system is them found through variational optimization. Once this is
% complete, properties of this ground state are plotted, such as the
% magnetization along z and the z-z correlation. Try changing the intensity
% of the field and the observables plotted.
%% Wipe Everything And Get Parameters
clear

%% Parameters
N = 30;		% Number of sites
D = 25;		% Bond dimension
g = 0.0;	% Transverse field strength
precision = 1e-8;

sigma = struct('x',[0,1;1,0],'y',[0,-1i;1i,0],'z',[1,0;0,-1],'plus',[0,1;0,0],'minus',[0,0;1,0],'id',eye(2));

%% Define Hamiltonian Operator
d = 2;
D_O = 5;
W = zeros(D_O,D_O,d,d);
W(1,1,:,:) = sigma.id;
W(2,1,:,:) = sigma.x;
W(3,1,:,:) = sigma.y;
W(4,1,:,:) = sigma.z;
W(5,1,:,:) = g*sigma.x;
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

%% Create Observables
identity = cell(1,N);
magnetization = cell(N,N);
correlation = cell(N,N);
for site = 1:N
    identity{site} = reshape(sigma.id,[1 1 d d]);
end

for site = 1:N
    magnetization(site,:) = identity;
    magnetization{site,site} = reshape(sigma.z,[1,1,d,d]);
    correlation{site,1} = reshape(sigma.z,[1 1 d d]);
    correlation{site,site} = reshape(sigma.z,[1 1 d d]);
end
correlation{1,1} = reshape(sigma.z*sigma.z,[1 1 d d]);

%% Do Optimization
state = sweep(randomMPS(N,D,d,1),{},-1);
[state,E,iter] = ground_search(state,H,100,precision,true);

%% Plot Observables
magn = real(expectationvalue(magnetization,state));
corr_zz = real(expectationvalue(correlation,state));

figure
plot(1:N,magn,'s--');
xlabel('$k$')
ylabel('$\langle\sigma^z_k\rangle$')
ylim([-1 1])
figure
hold on
plot(2:N,abs(corr_zz(2:end)))
xlabel('$k$')
ylabel('$|\langle\sigma^z_1\sigma^z_k\rangle|$')
set(gca,'yscale','log')
set(gca,'xscale','log')
