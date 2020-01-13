function settings = tdvp_settings(custom_settings)
% Define default settings
settings = struct;
settings.mode = 1;
settings.verbose = true;
settings.isreal = 0;
settings.customfun = [];
settings.integrator = struct;
settings.integrator.handle = @expmv;
settings.integrator.options = struct;
settings.integrator.options.tol = 1e-12;
settings.eigsolver = struct;
settings.eigsolver.handle = @eigs;
settings.eigsolver.mode = 'SR';
settings.eigsolver.options = struct;

% Merge settings with custom ones
if nargin == 1
    settings = mergestruct(settings,custom_settings);
end
% Update eigsolver settings for real matrices
if isequal(settings.eigsolver.handle,@eigs)
    settings.eigsolver.options.issym = true;
    settings.eigsolver.options.isreal = false;
    if settings.isreal
        settings.eigsolver.options.isreal = true;
        settings.eigsolver.mode = 'SA';
    end
end
end

function defaults = mergestruct(defaults,mods)
    labels = fieldnames(mods);
    for ind = 1:length(labels)
        l = labels{ind};
        if isstruct(mods.(l))
            defaults.(l) = mergestruct(defaults.(l),mods.(l));
        else
            defaults.(l) = mods.(l);
        end
    end
end
