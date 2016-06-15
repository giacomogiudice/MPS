% This script is very similar to 'density_mps' and compares the MPS 
% evolution to the exact dynamics of the Lindbladian. To be used only for
% small systems. Expected output is a plot of the error between the two
% methods, it should stabilize after a sufficient amount of time.
% Initialize all parameters with 'initialize'

%% More Parameters
D = d^2;
time = 0:dt:T;
t_sampling = 0:dt*n_sampling:T;
if mod(T,dt*n_sampling) % Add final time in case of incommensurability
    t_sampling(end+1) = T;
end

%% Initial MPS State
mps = cell(1,N);
vacuum = reshape(zeros(d),[1 1 D]);
vacuum(1,1,1) = 1;
for i = 1:N
    mps{i} = vacuum;
end

%% Construction of Evolution Operators
A = diag(sqrt(1:(d-1)),1);
I = eye(d);

% Identity 'state'
id = cell(1,N);
for i = 1:N
    id{i} = reshape(I,[1 1 D]);
end
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
direction = 1;
err = zeros(length(time)-1,1);
S = swap_matrix(d,N);
for step = 1:length(time)-1
    % Second order Trotter scheme
    mps = sweep(mps,U_odd_half,+direction,D_max,epsilon);
    mps = sweep(mps,U_even,-direction,D_max,epsilon);
    mps = sweep(mps,U_odd_half,+direction,D_max,epsilon);
    % Enforce Tr(rho) = 1
    if direction == +1
        mps{N} = mps{N}/braket(id,mps);
    else
        mps{1} = mps{1}/braket(id,mps);
    end
    direction = -direction;
    % Lindbladian evolution
    state = U_hat*state;
    err(step) = norm(S*expand_mps(mps) - state)/norm(state);
end
%% Plot Difference
figure
plot(time(2:end),err)
set(gca,'yscale','log')
xlabel('$t/\gamma$')
ylabel('$\| \rho_{\rm MPS} - \rho_{\rm LI} \|/\|\rho_{\rm LI}\|$')