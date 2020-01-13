function handle = fun_onesite(O,left,right)
% Creates a function handle to compute the application of the "effective
% two-site Hamiltonian". This is used internally for the 2-site TDVP
% routine. Diagrammatically this corresponds to the function
%
%                |  | |  |        | |
%    fun(M) = left -O-O- right = - M -
%                 \ | | /
%                    M
%
% Note: the index order for M is (left,first physical,second physical,right)
%
% INPUT
%   O:          Two MPO elements in ta cel array
%   left,right: partial contraction of the left and right blocks
%               respectively (see block functions)
% OUTPUT
%   handle:     function handle to compute the contraction for a given MPS
%               element

function W = new_element(M)
    W = ncon({left,M,O{1},O{2},right},{[-1,1,2],[1,3,5,7],[2,4,-2,3],[4,6,-3,5],[-4,7,6]});
end
handle = @new_element;
end
