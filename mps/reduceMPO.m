function U = reduceMPO(mpo_list)
% Combine different operators, and then reduce the bond dimension of the
% resulting operator by removing linear dependence within the tensors,
% similarly to arXiv:1611.02498.
%
% INPUT
%	mpo_list:	cell array of MPOs to combine (the order is important if
%				they do not commute)
% OUTPUT
%	U:			MPO corresponding to the multiplication of all the input 
%				operators

n_mpo = length(mpo_list);
N = length(mpo_list{1});
d = size(mpo_list{1}{1},3);
U = cell(1,N);
carryover = 1;
for site = 1:N
	% Combine all terms in list
	W = mpo_list{1}{site};
	for k = 2:n_mpo
		W = multiply(W,mpo_list{k}{site});
	end
	% Contract with previous carryover
	W = ncon({carryover,W},{[-1,1],[1,-2,-3,-4]});
	[chi_1,chi_2,~,~] = size(W);
	% Perform delinearization (first pass)
	M = reshape(permute(W,[1,3,4,2]),[chi_1*d^2,chi_2]);
	[R,colind] = rref(M);
	chi_new = length(colind);
	U{site} = permute(reshape(M(:,colind),[chi_1,d,d,chi_new]),[1,4,2,3]);
	carryover = R(1:chi_new,:); 
end

carryover = 1;
for site = N:(-1):1
	% Contract with previous carryover
	W = ncon({U{site},carryover},{[-1,1,-3,-4],[1,-2]});
	[chi_1,chi_2,~,~] = size(W);
	% Perform delinearization (second pass)
	M = reshape(permute(W,[2,3,4,1]),[chi_2*d^2,chi_1]);
	[R,colind] = rref(M);
	chi_new = length(colind);
	U{site} = permute(reshape(M(:,colind),[chi_2,d,d,chi_new]),[4,1,2,3]);
	carryover = R(1:chi_new,:).'; 
end
end