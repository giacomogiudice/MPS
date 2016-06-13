function psi_new = semidet(H,psi,dt,a,gamma,N)
% Semi-deterministic step used for Monte-Carlo integration. The
% deterministic part is integrated with 4th order Runge-Kutta, while the
% stochastic part is integrated with a simple Euler-Maruyama step
%
% INPUT
%   H:          matrix corresponding to the Hamiltonian
%   psi:        state vector to evolve
%   dt:         time step
%   a:          cell array of annihilation operators for each site
%   gamma:      coupling constant to the environment
%   N:          number of sites
% OUTPUT
%   psi_new:    state vector corresponding to the next step

% Runge Kutta order 4 on deterministic part
dt_half = 0.5*dt;
k_1 = -1i*H*psi;
k_2 = -1i*H*(psi + dt_half*k_1);
k_3 = -1i*H*(psi + dt_half*k_2);
k_4 = -1i*H*(psi + dt*k_3);
psi_new = psi + dt/6*(k_1 + 2*k_2 + 2*k_3 + k_4);

% Euler-Maruyama on stochastic part
dW = sqrt(gamma*dt)*randn(N,1);
for k = 1:N
    a_tilde = a{k}*psi;     % Pre-compute a_k|psi>
    dQ = 2*gamma*dt*real(psi'*a_tilde) + dW(k);
    psi_new = psi_new + dQ*a_tilde;
end

% normalize wave-vector
psi_new = psi_new/norm(psi_new);
end
