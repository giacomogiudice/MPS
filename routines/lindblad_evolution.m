% This script simulates the evolution of driven-dissipative Bose-Hubbard
% chains by computing the evolution of the density matrix at evenly spaced
% out time steps. in. The 'state' corresponds to the vectorized 
% density-matrix, and its evolution is determined by the non-unitary
% super-operator exp(L*t), where L is the Lindbladian.
% Initialize all parameters with 'initialize'

%% More Parameters
t_sampling = 0:dt*n_sampling:T;
if mod(T,dt*n_sampling) % Add final time in case of incommensurability
    t_sampling(end+1) = T;
end

%% Initialize State
rho = sparse(1);
vac = sparse(d,d);
vac(1,1) = 1;
for i = 1:N
    rho = kron(rho,vac);
end
state = reshape(rho,[],1);

%% Initialize Operators
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

U_hat = sparse(expm(L*(dt*n_sampling)));

%% Time Evolution
if showbar
    h = waitbar(0,'Computing...');
end
sample = 1;
% Compute observables at initial state
n = zeros(N,length(t_sampling));
g2 = zeros(N,N,length(t_sampling));
for i = 1:N
    n(i,sample) = real(trace(a{i}'*a{i}*rho));
end
for i = 1:N
    for j = 1:i
        g2(i,j,sample) = trace(a{j}'*a{i}'*a{j}*a{i}*rho)/(n(i)*n(j));
        g2(j,i,sample) = g2(i,j,sample);
    end
end

% Do evolution
for step = 1:length(t_sampling)-1
    state = U_hat*state;
    sample = sample + 1;
    rho = reshape(state,d^N,d^N);
    for i = 1:N
        n(i,sample) = real(trace(a{i}'*a{i}*rho));
    end
    for i = 1:N
        for j = 1:i
            g2(i,j,sample) = trace(a{j}'*a{i}'*a{j}*a{i}*rho)/(n(i,sample)*n(j,sample));
            g2(j,i,sample) = g2(i,j,sample);
        end
    end
    if showbar
        waitbar(step/length(t_sampling),h,sprintf('n = %.6f',max(real(n(:,sample)))));
    end
end
if showbar
    close(h)
end