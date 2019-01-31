% This script creates the MPO Hamiltonian corresponding to the Heisenberg
% model, with an external field of strength g along x. The ground state of
% this system is them found through variational optimization. Once this is
% complete, properties of this ground state are plotted, such as the
% magnetization along z and the z-z correlation. Try changing the intensity
% of the field and the observables plotted.
%% Wipe Everything And Get Parameters
clear

%% Parameters
N = 30;			% Number of sites
D = [8 16 24];	% Bond dimension
g = 0.22;		% Transverse field strength

sx = [0,1;1,0];
sy = [0,-1i;1i,0];
sz = [1,0;0,-1];
si = [1,0;0,1];

%% Define Hamiltonian Operator
d = 2;
D_O = 5;
W = zeros(D_O,D_O,d,d);
W(1,1,:,:) = si;
W(2,1,:,:) = sx;
W(3,1,:,:) = sy;
W(4,1,:,:) = sz;
W(5,1,:,:) = -g*sz;
W(5,2,:,:) = sx;
W(5,3,:,:) = sy;
W(5,4,:,:) = sz;
W(5,5,:,:) = si;

H = cell(1,N);
H{1} = W(D_O,:,:,:);
for i = 2:N-1
    H{i} = W;
end
H{N} = W(:,1,:,:);

%% Do Optimization
[state,output] = variational(H,D);

%% Plot Observables
magn = real(expectationvalue(state,sz));
corr_zz = real(expectationvalue(state,{sz,sz}));

figure(1)
plot(1:N,magn,'s--');
xlabel('$k$')
ylabel('$\langle\sigma^z_k\rangle$')
ylim([-1 1])

figure(2)
hold on
plot(1:(N-1),abs(corr_zz))
xlabel('$k$')
ylabel('$|\langle\sigma^z_1\sigma^z_k\rangle|$')
set(gca,'yscale','log')
set(gca,'xscale','log')
