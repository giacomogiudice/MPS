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
	W = ncon({left,M,O,right},{[-1,1,2],[1,4,3],[2,5,-3,3],[-2,4,5]});
end
handle = @new_element;
end
