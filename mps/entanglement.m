function [E,S] = entanglement(mps,direction)
% Computes the entanglement entropy between each partition of the MPS.

% INPUT
% 	mps:		cell array correponding to MPS representation of the state
%	direction:	direction the input MPS is canonized
% OUTPUT
%	E:			array corresponding to the entanglement at each partition
%	S:			cell array of the schmidt values at each partition
N = length(mps);
E = zeros(1,N-1);
S = cell(1,N-1);
% Go in the oppositie direction of the canonization
direction = -direction;
switch direction
	case +1 % Go left -> right
		M = mps{1};
		for site = 1:(N-1)
			[M,carryover,schmidt] = canonize(M,1,max(size(M,1),size(M,2)));
			S{site} = diag(schmidt);
			E(site) = -sum(S{site}.^2.*log(S{site}.^2));
			M = ncon({carryover,mps{site+1}},{[-1,1],[1,-2,-3]});
		end
	case -1 % Go right -> left
		M = mps{N};
		for site = N:(-1):2
			[M,carryover,schmidt] = canonize(M,-1,max(size(M,1),size(M,2)));
			S{site-1} = diag(schmidt);
			E(site-1) = -sum(S{site-1}.^2.*log(S{site-1}.^2));
			M = ncon({mps{site-1},carryover},{[-1,1,-3],[1,-2]});
		end
end
end

