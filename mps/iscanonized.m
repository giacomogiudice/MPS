function s = iscanonized(mps,direction,tolerance)
% Checks if an MPS is canonized along the direction specified. Optional 
% argument specifies the distance from the identity matrix which is
% acceptable for each site.
%
% INPUT
%   mps:		cell-array of rank-3 tensors representing the MPS 
%	direction:	specifies left (-1) or right (+1) canonization
%	tolerance:	(optional) acceptable error in canonical form. 
%				The canonical form must be no further distant (in  
%				Frobenius norm) to the identity than tol*D^2
% OUTPUT
%	s:			boolean which is set to true if condition is met

if nargin == 2
	tolerance = eps;
	
if ~iscell(mps)
	error('Expected cell array as first argument');
end

N = length(mps);
d = size(mps{1},3);

s = true;
blk = 1;
if direction == 1
	ind = 1:N;
elseif direction == -1
	ind = N:(-1):1;
else
	error('Unrecognized direction %d.',direction);
end

for site = ind
	blk = update_block(blk,mps{site},{},mps{site},direction);
	if norm(blk - eye(size(blk)),'fro') > tolerance*prod(size(blk))
		s = false;
		return
	end
end
end
