function [U,U_norm] = compressMPO(varargin)
% Combine different operators, and then compress the resulting operator
% by cutting the zeros. The norm is distributed evenly across all elements.
%
% INPUT
%	mps_1,mps_2,...:	variable number of MPOs to multiply together
%						(the order is important if they do not commute)
% OUTPUT
%	U:					MPO corresponding to the multiplication of all the 
%						input operators
%	U_norm:				Frobenius norm of the resulting operator

for n_mpo = 1:length(varargin)
	if ~iscell(varargin{n_mpo})
		break;
	end
end

N = length(varargin{1});
D = size(varargin{1}{1},3);
% Combine and reshape into an MPS
U = cell(1,N);
for site = 1:N
		W = varargin{1}{site};
	if n_mpo > 1
		for k = 2:n_mpo
			W = multiply(W,varargin{k}{site});
		end
	end
	[b_1,b_2,~,~] = size(W);
	U{site} = reshape(W,[b_1,b_2,D*D]);
end
% Compress by cutting down to machine precision
[U,U_norm] = sweep(U,{},-1);
[U,U_sign] = sweep(U,{},1,1e6,eps);
U_norm = real(U_norm*U_sign);
% Distribute the norm along all sites
for site = 1:N
	U{site} = U_norm^(1/N)*U{site};
end
% Reshape each element to a rank-4 tensor
for site = 1:N
	[b_1,b_2,~,~] = size(U{site});
	U{site} = reshape(U{site},[b_1,b_2,D,D]);
end
end