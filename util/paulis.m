function [sx,sy,sz,si] = paulis()
% PAULIS Returns the Pauli matrices
%
% [sx,sy,sz,si] = paulis()
%

sx = [0,1;1,0];
sy = [0,-1i;1i,0];
sz = [1,0;0,-1];
si = [1,0;0,1];
end
