function mpo = reduceMPO(mpo_list,iter_max)
% Combine different operators, and then reduce the bond dimension of the
% resulting operator by removing linear dependence within the tensors,
% similarly to arXiv:1611.02498.
%
% INPUT
%	mpo_list:	cell array of MPOs to combine (the order is important if
%				they do not commute)
% OUTPUT
%	mpo:			MPO corresponding to the multiplication of all the input 
%				operators
if nargin < 2
	iter_max = 10;
end
n_mpo = length(mpo_list);
N = length(mpo_list{1});
d = size(mpo_list{1}{1},3);
mpo = cell(1,N);


for iter = 1:iter_max
	carryover = 1;
	improvement = 0;
	for site = 1:N
		% Combine all terms in list
		if iter == 1
			W = mpo_list{1}{site};
			for k = 2:n_mpo
				W = multiply(W,mpo_list{k}{site});
			end
		else
			W = mpo{site};
		end
		% Contract with previous carryover
		W = ncon({carryover,W},{[-1,1],[1,-2,-3,-4]});
		[chi_1,chi_2,~,~] = size(W);
		% Perform delinearization (first pass)
		M = reshape(permute(W,[1,3,4,2]),[chi_1*d^2,chi_2]);
		[R,colind] = rref(M);
		chi_new = length(colind);
		mpo{site} = permute(reshape(M(:,colind),[chi_1,d,d,chi_new]),[1,4,2,3]);
		carryover = R(1:chi_new,:); 
		improvement = improvement - chi_new + chi_2;
	end
	carryover = 1;
	for site = N:(-1):1
		% Contract with previous carryover
		W = ncon({mpo{site},carryover},{[-1,1,-3,-4],[1,-2]});
		[chi_1,chi_2,~,~] = size(W);
		% Perform delinearization (second pass)
		M = reshape(permute(W,[2,3,4,1]),[chi_2*d^2,chi_1]);
		[R,colind] = rref(M);
		chi_new = length(colind);
		mpo{site} = permute(reshape(M(:,colind),[chi_2,d,d,chi_new]),[4,1,2,3]);
		carryover = R(1:chi_new,:).'; 
		% [chi_1,chi_new]
		improvement = improvement - chi_new + chi_1;
	end
	if improvement == 0
		break;
	end
end
end