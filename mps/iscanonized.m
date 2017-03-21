function s = iscanonized(mps,direction,varargin)
% Checks if an MPS is canonized along the direction specified. Optional 
% argument specifies the distance from the identity matrix which is
% acceptable for each site.
%
% INPUT
%   mps:        cell-array of rank-3 tensors representing the MPS 
%   direction:  specifies left (-1) or right (+1) canonization
%   tolerance:  (optional) acceptable error in distance from identity
% OUTPUT
%   s:          boolean which is set to true if condition is met

tolerance = 1e-10;
if ~isempty(varargin)
    tolerance = varargin{1};
end
    
if ~iscell(mps)
    disp('Expected cell array as first argument');
    return
end

N = length(mps);
d = size(mps{1},3);

if direction == 1
    for i = 1:N
        s = 0;
        for sigma= 1:d
            s = s + mps{i}(:,:,sigma)'*mps{i}(:,:,sigma);
        end
        if norm(s - eye(length(s))) > tolerance
            s = false;
            return
        end
    end
elseif direction == -1
    for i = 1:N
        s = 0;
        for sigma= 1:d
            s = s + mps{i}(:,:,sigma)*mps{i}(:,:,sigma)';
        end
        if norm(s - eye(length(s))) > tolerance
            s = false;
            return
        end
    end   
end
s = true;
end
