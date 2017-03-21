% This stript tests the first order Trotter decomposition

close all
clear all

%% Initial Parameters
N = 8;
J = 1;
U = 1.5;
d = 2;
dt = 0.1; % Time increment
test_tolerance = 1e-12;

%% Two-site Operator

Ham_pair = randn(d^2,d^2);
Ham_pair = Ham_pair'*Ham_pair;

%% 1st Order Trotter Splitting
U_odd = cell(1,N);
U_even = cell(1,N);
[V,W] = trotter(expm(-1i*dt*Ham_pair),d);
for site = 1:2:N-1
	U_odd{site} = V;
	U_odd{site+1} = W;
end
for site = 2:2:N-1
	U_even{site} = V;
	U_even{site+1} = W;
end

U_even{1} = reshape(eye(2),[1 1 d d]);
if mod(N,2)
	U_odd{N} = reshape(eye(2),[1 1 d d]);
else
	U_even{N} = reshape(eye(2),[1 1 d d]);
end

%% Kron Version
 H_odd_kron = 0;

 for i = 1:2:N-1
    H_odd_kron = H_odd_kron + kron(kron(eye(d^(i-1)),Ham_pair),eye(d^(N-i-1)));
 end
 U_odd_kron = expm(-1i*dt*H_odd_kron);

 H_even_kron = 0;
for i = 2:2:N-1
    H_even_kron = H_even_kron + kron(kron(eye(d^(i-1)),Ham_pair),eye(d^(N-i-1)));
end
U_even_kron = expm(-1i*dt*H_even_kron);

%% Compare
fprintf('Testing U_odd...');
assert(norm(expandMPO(U_odd) - U_odd_kron) < test_tolerance);
fprintf('\tdone\n');
fprintf('Testing U_even...');
assert(norm(expandMPO(U_even) - U_even_kron) < test_tolerance);
fprintf('\tdone\n');

fprintf('All tests passed!\n');
