function handle = fun_zerosite(left,right)
% Creates a function handle to compute the application of the "effective
% zero-site Hamiltonian". This is used internally for the TDVP routine to
% evolve the carryover backward in time. Diagrammatically this corresponds
% to the function
%                  
%                |   |
%    fun(C) = left - right = -W-
%                 \ /
%                  C	
%
% INPUT
%	left,right:	partial contraction of the left and right blocks 
%				respectively (see block functions)
% OUTPUT
%	handle:		function handle to compute the contraction for a given 
%				carryover

function W = new_element(C)
	W = ncon({left,C,right},[-1,1,3],[1,2],[-2,2,3]);
end
handle = @new_element;
end
