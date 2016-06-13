% An alternative time evolution that distributes data along the GPU.
% See 'density_mps' for further information.

%% Function Definitions
gpuComplex = @(x) gpuArray(complex(x));
gpuCell = @(x) cellfun(gpuComplex,x,'UniformOutput',false);
gpuCellGather = @(x) cellfun(@gather,x,'UniformOutput',false);

%% Load on GPU
U_odd = gpuCell(U_odd);
U_even = gpuCell(U_even);
U_odd_half = gpuCell(U_odd_half);

%% Time Evolution
if showbar, h = waitbar(0,'Sweeping...'); end
sample = 1;
n = zeros(N,length(t_sampling));
g2 = zeros(N,N,length(t_sampling));
bond_size = zeros(1,length(t_sampling));
% Calculate observables at initial step
rho = gpuCellGather(state);
n(:,sample) = occupation_mps(occupation,rho);
g2(:,:,sample) = correlation_mps(correlation,rho,n(:,sample));
bond_size(sample) = maxbond(rho);
direction = -1;
% Initial sweep
state = sweep(state,U_odd_half,+direction,D_max,epsilon);

for step = 1:length(time)-2
    % Second order Trotter scheme
    if mod(step,n_sampling) == 0 % Time to measure
        sample = sample + 1;
        state = sweep(state,U_even,-direction,D_max,epsilon);
        state = sweep(state,U_odd_half,+direction,D_max,epsilon);
        rho = gpuCellGather(state);
        % Enforce Tr(rho) = 1
        if direction == +1
            rho{N} = rho{N}/braket(id,rho);
        else
            rho{1} = rho{1}/braket(id,rho);
        end
        % Calculate observables
        n(:,sample) = occupation_mps(occupation,rho);
        g2(:,:,sample) = correlation_mps(correlation,rho,n(:,sample));
        bond_size(sample) = maxbond(rho);
        % Do last sweep to return to initial position
        state = sweep(state,U_odd_half,-direction,D_max,epsilon);
        direction = -direction;
        if showbar
            waitbar(step/length(time),h,sprintf('n = %.6f',max(real(n(:,sample)))));
        end
    else % No measurements, go faster
        state = sweep_mult(state,{U_even U_odd},-direction,D_max,epsilon);
        direction = -direction;
    end
end
%  Last step
sample = sample + 1;
state = sweep(state,U_even,-direction,D_max,epsilon);
state = sweep(state,U_odd_half,+direction,D_max,epsilon);
% Enforce Tr(rho) = 1
rho = gpuCellGather(state);
if direction == +1
    rho{N} = rho{N}/braket(id,rho);
else
    rho{1} = rho{1}/braket(id,rho);
end
% Calculate observables
n(:,sample) = gather(occupation_mps(occupation,rho));
g2(:,:,sample) = gather(correlation_mps(correlation,rho,n(:,sample)));
bond_size(sample) = maxbond(rho);
if showbar
    close(h);
end

%% Gather All Distributed Arrays
state = rho;
U_odd = gpuCellGather(U_odd);
U_even = gpuCellGather(U_even);
U_odd_half = gpuCellGather(U_odd_half);