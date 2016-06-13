% This script uses the MPS method to simulate the dynamics of the density
% matrix in a driven-dissipative Bose-Hubbard chain. The 'state'
% corresponds to the vectorized density-matrix, and its evolution is
% determined by the non-unitary super-operator exp(L*t), where L is the
% Lindbladian. This is then handled in MPO form.
% Initialize all parameters with 'initialize'

%% More Parameters
D = d^2;
time = 0:dt:T;
t_sampling = 0:dt*n_sampling:T;
if mod(T,dt*n_sampling) % Add final time in case of incommensurability
    t_sampling(end+1) = T;
end

%% Initial MPS State
if ~iscell(state)
    state = cell(1,N);
    vacuum = reshape(zeros(d),[1 1 D]);
    vacuum(1,1,1) = 1;
    for i = 1:N
        state{i} = vacuum;
    end
end

%% Construction of Evolution Operators
A = diag(sqrt(1:(d-1)),1);
I = eye(d);
% 2-site interaction Lindbladian
L_int = 1i*J*( ...
    + kron(kron(I,A),kron(I,A')) + kron(kron(I,A'),kron(I,A)) ...
    - kron(kron(A',I),kron(A,I)) - kron(kron(A,I),kron(A',I)) ...
    );
% 1-site Hamiltonian
H_onsite = d_omega*(A'*A) + U*A'*(A'*A)*A + F*(A' + A);
% 1-site Lindbladian
L_onsite = 1i*(kron(H_onsite',I) - kron(I,H_onsite)) ...
    - gamma/2*(kron(I,A'*A) + kron(A'*A,I) - 2*kron(A,A));
% 2-site Lindbladian of onsite interaction and dissipation
L_twosite = 1i*(kron(kron(H_onsite,I),kron(I,I)) + kron(kron(I,I),kron(H_onsite,I)) ...
    - kron(kron(I,H_onsite),kron(I,I)) - kron(kron(I,I),kron(I,H_onsite))) ...
    - gamma/2*(kron(kron(I,I),kron(A'*A,I)) + kron(kron(I,I),kron(I,A'*A)) + ...
    + kron(kron(A'*A,I),kron(I,I)) + kron(kron(I,A'*A),kron(I,I)) ...
    - 2*kron(kron(A,A),kron(I,I)) - 2*kron(kron(I,I),kron(A,A)));
% Decompose Lindbladian as MPOs
[U_odd_half,~] = trotter(expm((L_int + L_twosite/2)*(dt/2)),N,D);
[U_odd,U_even] = trotter(expm((L_int + L_twosite/2)*dt),N,D);

% Fix boundaries
U_even_half{1} = reshape(expm(L_onsite*dt/4),[1 1 D D]);
U_even{1} = reshape(expm(L_onsite*dt/2),[1 1 D D]);
if mod(N,2)
    U_odd_half{N} = reshape(expm(L_onsite*dt/4),[1 1 D D]);
    U_odd{N} = reshape(expm(L_onsite*dt/2),[1 1 D D]);
else
%     U_even_half{N} = reshape(expm(L_onsite*dt/4),[1 1 D D]);
    U_even{N} = reshape(expm(L_onsite*dt/2),[1 1 D D]);
end

%% Measurement Operators
id = cell(1,N);
occupation = cell(1,N);
correlation = cell(N,N);
for i = 1:N
    id{i} = reshape(I,[1 1 D]);
end
for i = 1:N
    occupation{i} = id;
    occupation{i}{i} = reshape(A'*A,[1 1 D]);
end
for i = 1:N
    for j = 1:i
        correlation{i,j} = id;
        if i == j
            correlation{i,j}{i} = reshape(A'*(A'*A)*A,[1 1 D]);
        else
            correlation{i,j}{i} = reshape(A'*A,[1 1 D]);
            correlation{i,j}{j} = reshape(A'*A,[1 1 D]);
        end
    end
end
maxbond = @(s) max(cellfun(@(x) max(size(x,1),size(x,2)),s));

%% Switch to GPU If Available
if gpuDeviceCount
    dm_gpu;
    return
end

%% Time Evolution
if showbar, h = waitbar(0,'Sweeping...'); end
sample = 1;
n = zeros(N,length(t_sampling));
g2 = zeros(N,N,length(t_sampling));
bond_size = zeros(1,length(t_sampling));
% Calculate observables at initial step
n(:,sample) = occupation_mps(occupation,state);
g2(:,:,sample) = correlation_mps(correlation,state,n(:,sample));
bond_size(sample) = maxbond(state);
direction = -1;
% Initial sweep
state = sweep(state,U_odd_half,+direction,D_max,epsilon);

for step = 1:length(time)-2
    % Second order Trotter scheme
    if mod(step,n_sampling) == 0 % Time to measure
        sample = sample + 1;
        state = sweep(state,U_even,-direction,D_max,epsilon);
        state = sweep(state,U_odd_half,+direction,D_max,epsilon);
        % Enforce Tr(rho) = 1
        if direction == +1
            state{N} = state{N}/braket(id,state);
        else
            state{1} = state{1}/braket(id,state);
        end
        % Calculate observables
        n(:,sample) = occupation_mps(occupation,state);
        g2(:,:,sample) = correlation_mps(correlation,state,n(:,sample));
        bond_size(sample) = maxbond(state);
        % Do last sweep to return to initial position
        state = sweep(state,U_odd_half,-direction,D_max,epsilon);
        direction = -direction;
        if showbar
            waitbar(step/length(time),h,sprintf('n = %.6f',max(real(n(:,sample)))));
        end
    else % No measurements, go faster
        state = sweep(state,U_even,-direction,D_max,epsilon);
        state = sweep(state,U_odd,+direction,D_max,epsilon);
    end
end
%  Last step
sample = sample + 1;
state = sweep(state,U_even,-direction,D_max,epsilon);
% state = sweep(state,U_odd_half,+direction,D_max,epsilon);
% Enforce Tr(rho) = 1
if direction == +1
    state{N} = state{N}/braket(id,state);
else
    state{1} = state{1}/braket(id,state);
end
% Calculate observables
n(:,sample) = occupation_mps(occupation,state);
g2(:,:,sample) = correlation_mps(correlation,state,n(:,sample));
bond_size(sample) = maxbond(state);
if showbar
    close(h);
end
