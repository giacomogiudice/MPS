% This script computes the steady-state of driven-dissipative Bose-Hubbard
% chains by inverting the Lindbladian.
% Initialize all parameters with 'initialize'

%% Initialize State and Operators
vac =  reshape(sparse(ground(N,d)*ground(N,d)'),[],1);
a = annihilation(N,d);

H = sparse(d^N,d^N);
for i = 1:N
    % omega, U, F
    H = H + d_omega*(a{i}'*a{i})...
        + U*a{i}'*a{i}'*a{i}*a{i}...
        + F*(a{i}' + a{i});
    % J term
    if i ~= N
        H = H - J*a{i}*a{i+1}' - J*a{i}'*a{i+1};
    end
end

I = speye(d^N);
L = 1i*(kron(H',I) - kron(I,H));
for i = 1:N
    L = L - 0.5*gamma*(kron(I,a{i}'*a{i}) + kron(a{i}'*a{i},I) - 2*kron(a{i},a{i}));
end

%% Compute Inversion
L = L + vac*reshape(I,1,[]);
rho = reshape(L\vac,d^N,d^N);

%% Compute Observables
n = zeros(N,1);
g2 = zeros(N,N);
for i = 1:N
    n(i) = real(trace(a{i}'*a{i}*rho));
end
for i = 1:N
    for j = 1:i
        g2(i,j) = trace(a{j}'*a{i}'*a{j}*a{i}*rho)/(n(i)*n(j));
        g2(j,i) = g2(i,j);
    end
end
