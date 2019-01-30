function settings = variational_settings(custom_settings)
% This function creates a struct corresponding to the settings that are
% passed to 'variational.m'. Optional settings are passed to this function
% and it completes the struct with the default settings. 
% The list of valid fields one can pass include:
%	maxit 			Maximum numer of iterations
% 	tol:			halting tolerance, i.e. when 
%					|energy[k] - energy[k-1]| < precision*max(|energy[k]|,1) 
% 	isreal:			constrain the tensors to be real (not tested)
% 	verbose: 		print output at each iteration
% 	initial: 		structure containing the initial guess. An MPS should 
%					be provided in settings.initial.mps
%					(WARNING: all elements must be left (+1) canonized!) 
% 	orthogonalize: 	cell array of MPS that should be orthogonal to the optimized MPS.
%					(WARNING: all elements must be left canonized!)  
% 	eigsolver:		structure of options for the eigensolver. This includes:
%					eigsolver.handle:	function pointer to eigensolver ('eigs' by default)
%					eigsolver.mode:		search for smallest or largest eigenvalue
%					eigsolver.options:	internal options of eigensolver
%
% INPUT
%	custom_settings:	struct corresponding to the options for the routine
% OUTPUT
%	settings:			full struct of settings, including defaults

% Define default settings
settings = struct;
settings.maxit = 10;
settings.tol = 1e-9;
settings.isreal = false;
settings.verbose = true;
settings.initial = struct;
settings.orthogonalize = {};
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
% Do some minimal checks
assert(mod(settings.maxit,1) == 0,'%s must be an integer.','settings.maxit');
assert(settings.tol <= 1 & settings.tol >= 0,'%s must be between 0 and 1.','settings.tol');

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
