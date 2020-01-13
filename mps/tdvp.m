function varargout = tdvp(mps,mpo,timesteps,settings)

% [mps,output,stats] = tdvp(mps,mpo,timesteps,settings)
if ~exist('settings','var') || isempty(settings)
    settings = tdvp_settings();
else
    settings = tdvp_settings(settings);
end

N = length(mps);
niter = length(timesteps);

% Initialize block storage
blocks = cell(1,N+1);
blocks{1} = 1;
blocks{N+1} = 1;
for site = N:(-1):2
    blocks{site} = update_block(blocks{site+1},mps{site},mpo{site},mps{site},-1);
end
energy_prev = update_block(blocks{2},mps{1},mpo{1},mps{1},-1);

if settings.mode == 1
    take_step = @tdvp_step;
elseif settings.mode == 2
    take_step = @tdvp2_step;
else
    error('Unavailable mode');
end

savestats = false;
customeval = false;
if nargout > 2
    savestats = true;
    stats.energy = zeros(1,niter);
    stats.energydiff = zeros(1,niter);
    if isfield(settings,'customfun') && isa(settings.customfun,'function_handle')
        customeval = true;
        stats.customvalue = cell(1,niter);
    end
end

if settings.verbose
    fprintf('Iter\t      Energy\t Energy Diff\tLap Time [s]\tBond Dim\n')
end
for iter = 1:niter
    if settings.verbose, tic; end
    [mps,blocks,energy] = take_step(mps,mpo,timesteps(iter),blocks,settings);
    laptime = toc;
    energydiff = energy - energy_prev;
    maxbond = max(cellfun(@(t) max(size(t,1),size(t,2)),mps));
    % Print results of iteration
    if settings.verbose
        fprintf('%4d\t%12g\t%12g\t%12.1f\t%8d\n',iter,energy,energydiff,laptime,maxbond);
    end
    % Compute observables
    if savestats
        stats.energy(iter) = energy;
        stats.energydiff(iter) = energydiff;
        if customeval
            try
                stats.customvalue{iter} = settings.customfun(mps,mpo,settings);
            catch me
                warning('Failed to evaluate custom function.');
                disp(getReport(me,'extended','hyperlinks','on'));
                stats.customfun{iter} = [];
            end
        end
    end
    energy_prev = energy;
end
output.iter = iter;
output.energy = energy;

varargout{1} = mps;
varargout{2} = output;
if savestats
    varargout{3} = stats;
end
end
