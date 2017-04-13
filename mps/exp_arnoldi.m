function [v_out,mk] = exp_arnoldi(v_in,fun,dt,m,tolerance)
% This function computes an approximation exp(A*dt)*v_in, for a matrix A.
% It does not compute the full exponential, but instead uses a 
% Krylov-space technique that only requires the application of A to a 
% vector. However, A is not provided as a matrix, but a function 
% fun(v) = A*v must be provided instead. Additionally, v_in does not have
% to be a vector, but can be a tensor of any rank, as long as fun(v_in) is
% a tensor of the same rank. This is an adapted version of 'expv' from the
% Expokit package.
%
% INPUT
%	v_in:		input tensor (or pseudo-vector) corresponding to the state
%				one wants to evolve
%   fun:		function that computes the application of a pseudo-matrix 
%				to a pseudo-vector
%   dt:			time step (can be imaginary)
%	m:			(optional) maximum number of additional Krylov 
%				pseudo-vectors to be stored
%	tolerance:	(optional) stopping condition for the number of Krylov 
%				vectors to store, measured as the the norm of the next
%				Krylov vector (see Expokit paper by R. B. Sidje)
% OUTPUT
%	v_out:		output tensor (or pseudo-vector) corresponding to the 
%				evolved state
%	mk:			number of additional Krylov pseudo-vectors effectively used

% Handle optional arguments
if nargin < 5
	tolerance = 1e-6;
end
if nargin < 4
	m = min(20,numel(v_in));
end

scalar_product = @(x,y) conj(reshape(x,1,[]))*reshape(y,[],1);
mred = 0;

hessenberg = zeros(m+2,m+2);
krylov = cell(1,m+1);
normv = sqrt(scalar_product(v_in,v_in));
krylov{1} = v_in/normv;

% Build Krylov subspace and diagonalize with Arnoldi
for j = 1:m
	work = fun(krylov{j});
	for k = 1:j
		hessenberg(k,j) = scalar_product(krylov{k},work);
		work = work - hessenberg(k,j)*krylov{k};
	end
	hessenberg(j+1,j) = sqrt(scalar_product(work,work));
	if hessenberg(j+1,j) < tolerance
		% Happy ending, the precision is sufficient 
		mred = j;
		break
	end
	krylov{j+1} = work/hessenberg(j+1,j);
end
if mred == 0
	hessenberg(m+2,m+1) = 1;
	mred = m + 2;
end

% Compute exponential approximation in Krylov subspace
U = expm(dt*hessenberg(1:mred,1:mred));
U = normv*U(:,1);

% Build output vector
mk = min(mred,m+1);
v_out = U(1)*krylov{1};
for j = 2:mk
	v_out = v_out + U(j)*krylov{j};
end
end