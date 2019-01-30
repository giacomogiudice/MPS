function [M,energy] = optimization_step(M,fun,settings)
% Options for eigs routine
if nargin >= 3
	eigsolver = settings.eigsolver.handle;
	eigsolver_options = settings.eigsolver.options;
	eigsolver_mode = settings.eigsolver.mode;
else
	eigsolver = @eigs;
	eigsolver_options.issym = true;
	eigsolver_options.isreal = false;
	eigsolver_mode = 'sr';
end
eigsolver_options.v0 = reshape(M,[],1);
d_M = size(M);
v_fun = @(v) reshape(fun(reshape(v,d_M)),[],1);
% Find smallest real eigenvalue
try
	[M,energy] = eigsolver(v_fun,prod(d_M),1,eigsolver_mode,eigsolver_options);
catch
	warning('mps:variational:optimizationfail','Cannot find lowest eigenvalue. Moving on.');
	M_prime = fun(M);
	energy = M_prime(:)'*M(:);
end
M = reshape(M,d_M);
end