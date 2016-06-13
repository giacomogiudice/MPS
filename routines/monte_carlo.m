% This script simulates the evolution of driven-dissipative Bose-Hubbard
% chains by integrating the Stochastic Schrodinger equation. Thus the
% outputted state will be a realization of this process.
% Initialize all parameters with 'initialize'

%% Initialize State
scheme = @semidet;

state = ground(N,d);
time = 0:dt:T;
t_sampling = 0:dt*n_sampling:T;
if mod(T,dt*n_sampling) % Add final time in case of incommensurability
    t_sampling(end+1) = T;
end

%% Initialize Operators
a = annihilation(N,d);
occupation = cell(1,N);
correlation = cell(N,N);

for i = 1:N
    occupation{i} = a{i}'*a{i};
    for j = 1:i-1
        correlation{i,j} = a{j}'*a{i}'*a{j}*a{i};
    end
    correlation{i,i} = a{i}'*a{i}'*a{i}*a{i};
end

%% Initialize Hamiltonian
H = sparse(d^N,d^N);
for i = 1:N
    % omega, U, F and non-unitary term
    H = H + (d_omega - 0.5*1i*gamma)*(a{i}'*a{i})...
        + U*a{i}'*a{i}'*a{i}*a{i}...
        + F*(a{i}' + a{i});
    % J term
    if i ~= N
        H = H - J*a{i}*a{i+1}' - J*a{i}'*a{i+1};
    end
end

%% Time Evolution
if showbar
    h = waitbar(0,'Computing...');
end
sample = 1;
n = zeros(N,length(t_sampling));
g2 = zeros(N,N,length(t_sampling));

% Compute observables at initial state
n(:,sample) = occupation_mc(occupation,state);
g2(:,:,sample) = correlation_mc(correlation,state,n(:,sample));
% Do evolution
for step = 1:length(time)-1
    state = scheme(H,state,dt,a,gamma,N);
    if mod(step,n_sampling) == 0 || step == length(time)
        sample = sample + 1;
        n(:,sample) = occupation_mc(occupation,state);
        g2(:,:,sample) = correlation_mc(correlation,state,n(:,sample));
        if showbar
            waitbar(step/length(time),h,sprintf('n = %.6f',max(real(n(:,sample)))));
        end
    end
end
if showbar
    close(h)
end