function [M,energy] = optimization_step(M,fun)
	d_M = size(M);
	v_fun = @(v) reshape(fun(reshape(v,d_M)),[],1);
	% Options for eigs routine
	options.isreal = 0;
	options.issym = 1;
	options.v0 = reshape(M,[],1);
	% Find smallest real eigenvalue
	try
		[M,energy] = eigs(v_fun,prod(d_M),1,'sr',options);
	catch
		try
			options.v0 = [];
			[M,energy] = eigs(v_fun,prod(d_M),1,'sr',options);
		catch
			warning('mps:optimizationfail','Cannot find lowest eigenvalue. Moving on.');
			M_prime = fun(M);
			energy = M_prime(:)'*M(:);
		end
	end
	M = reshape(M,d_M);
end