function [C,carryover] = canonize_fast(M,direction)
% Computes canonization of MPS site M according
% to left- (direction == +1) or right- canonization (direction == -1).
% This is faster than 'canonize' since it uses QR decomposition instead of
% SVD.
%
% INPUT
%	M:			site element of MPS, represented as a rank-3 tensor
%	direction:	specifies left (-1) or right (+1) canonization\
% OUTPUT
%	C:			MPS element canonized, A or B
%	carryover:	product to multiply with next neighbor, S*V' or U*S

[b_1,b_2,d] = size(M);
switch direction
	case +1 % Going right
		C = reshape(permute(M,[1,3,2]),d*b_1,b_2);
		[Q,R] = qr(C,0);
		C = permute(reshape(Q,[b_1,d,size(Q,2)]),[1,3,2]);
		carryover = R;
		
	case -1 % Going left
		C = reshape(M,[b_1,b_2*d]);
		[Q,R] = qr(C',0);
		C = reshape(Q',[size(Q,2),b_2,d]);
		carryover = R';
end
end
