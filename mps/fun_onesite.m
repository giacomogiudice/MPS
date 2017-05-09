function handle = fun_onesite(O,left,right)
% Creates a function handle to compute the application of the "effective
% one-site Hamiltonian". This is used internally for the TDVP and ground 
% search routines. Diagrammatically this corresponds to the function
%                          
%                |  |  |        |
%    fun(M) = left -O- right = -W-
%                 \ | /
%                   M	
%
% INPUT
%	O:			MPO element in the contraction	
%	left,right:	partial contraction of the left and right blocks 
%				respectively (see block functions)
% OUTPUT
%	handle:		function handle to compute the contraction for a given MPS
%				element

function W = new_element(M)
	W = contract(right,3,2,M,3,2);
	W = contract(W,4,[2,4],O,4,[2,4]);
	W = contract(left,3,[2,3],W,4,[2,3]);
end
handle = @new_element;
end
