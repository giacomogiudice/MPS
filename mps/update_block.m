function new_block = update_block(prev_block,M_1,O,M_2,direction)
% To be used in the block abstraction when contracting an MPS. This function
% computes the update of the new block, in the direction specified.
% All blocks are rank-3 or rank-2 (depending if the operator has a virtual 
% bond dimension) and the indices are ordererd with the convention 
% (bra,ket,operator). The corresponding tensor diagrams are
%
%   -M_1                                M_1-
%    |  \                              /  |
%   -O-- prev_block    or    prev_block --O-
%    |  /                              \  |
%   -M_2                                M_2-
%   
% direction == -1(right) or   direction == +1 (left)
%
% INPUT
%	prev_block:	rank-3 tensor corresponding to the partial contraction
%				of all the blocks to the right or the left
%   M_1:      	MPS element corresponding to the bra-state
%	O:			MPO element to sandwich between the MPS elements, an empty 
%				element can be provided
%   M_2:      	MPS element corresponding to the ket-state
%	direction:	direction in which the block are being updated, specify
%				left (-1) or right (+1)
% OUTPUT
% 	new_block:	rank-3 tensor corresponding to the previous block contracted
%				with the next layer

switch direction
	case +1 % left -> right
		if isempty(O)
			new_block = ncon({prev_block,conj(M_1),M_2},{[1,2],[1,-1,3],[2,-2,3]});
		else
			new_block = ncon({prev_block,conj(M_1),O,M_2},{[1,2,3],[1,-1,4],[3,-3,4,5],[2,-2,5]});
		end
	case -1 % right -> left
		if isempty(O)
			new_block = ncon({prev_block,conj(M_1),M_2},{[1,2],[-1,1,3],[-2,2,3]});
		else
			new_block = ncon({prev_block,conj(M_1),O,M_2},{[1,2,3],[-1,1,4],[-3,3,4,5],[-2,2,5]});
		end
end
end
