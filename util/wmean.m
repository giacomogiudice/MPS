function y = wmean(x,w,dim)
% Takes weighted average or mean value of x using non-negative weights w 
% along the dimension dim. Each element of X requires a corresponding
% weight, and hence the size of w must match that of x.
% 
% INPUT
%   x:		array or matrix containing the values to average
%	w:		weights associated to each element x
%	dim:	dimension on which to take the mean value
% OUTPUT
%   y:      column vector containing the weighted mean

if nargin<2
    error('Not enough input arguments.');
end

% Check that dimensions of X match those of W.
if(~isequal(size(x), size(w)))
    error('Inputs x and w must be the same size.');
end

% Check that all of W are non-negative.
if (any(w(:)<0))
    error('All weights, W, must be non-negative.');
end

% Check that there is at least one non-zero weight.
if (all(w(:)==0))
    error('At least one weight must be non-zero.');
end

if nargin==2, 
  % Determine which dimension SUM will use
  dim = min(find(size(x)~=1));
  if isempty(dim), dim = 1; end
end

y = sum(w.*x,dim)./sum(w,dim);