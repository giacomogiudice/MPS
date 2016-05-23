function s = iscanonized(mps,direction,varargin)

tolerance = 1e-4;
if ~isempty(varargin)
    tolerance = varargin{1};
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