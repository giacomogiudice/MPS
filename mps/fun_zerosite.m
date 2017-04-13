function handle = fun_zerosite(left,right)
% Creates a function handle to compute the application of the "effective
% zero-site Hamiltonian". This is used internally for the TDVP routine to
% evolve the carryover backward in time. Diagrammatically this corresponds
% to the function
%                 -C-
%                /   \
%   fun(C) = left --- right
%                \   /
%                 -C-	
%
% INPUT
%	left,right:	partial contraction of the left and right blocks 
%				respectively (see block functions)
% OUTPUT
%	handle:		function handle to compute the contraction for a given 
%				carryover

function W = new_element(C)
	W = contract(right,3,2,C,2,2);
	W = contract(left,3,[2,3],W,3,[3,2]);
end
handle = @new_element;
end
