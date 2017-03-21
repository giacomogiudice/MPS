function [C,carryover] = canonize(M,direction,varargin)
% Computes canonization of MPS site M according
% to left- (direction == +1) or right- canonization (direction == -1).
% It then optionally performs a decimation if two additional are specified
% The decimation is computed by taking the sum of the square of the singular
% values and then keeping only the ones under (1 - epsilon)
%
% INPUT
%   M:          site element of MPS, represented as a rank-3 tensor
%   direction:  specifies left (-1) or right (+1) canonization
%   D_max:      (optional) maximum allowed bond length
%   epsilon:    (optional) tolerated error in decimation
% OUTPUT
%   C:          MPS element canonized, A or B
%   carryover:  product to multiply with next neighbor, S*V' or U*S

if isempty(varargin)
    [C,carryover] = canonize_fast(M,direction);
    return
else
    D_max = varargin{1};
    epsilon = varargin{2};
end

[b_1,b_2,d] = size(M);

switch direction
    case +1 % Going right
        C = reshape(permute(M, [1,3,2]),d*b_1,b_2);
        [U,S,V] = svd(C,'econ');
        if ~isempty(D_max)
            [U,S,V,D_cut] = trim(U,S,V,D_max,epsilon);
            % reshape V'
            C = permute(reshape(U,[b_1,d,D_cut]),[1,3,2]);
        else
            C = permute(reshape(U,[b_1,d,size(U,2)]),[1,3,2]);
        end
        carryover = S*V';
        
    case -1 % Going left
        C = reshape(M,[b_1,b_2*d]);
        [U,S,V] = svd(C,'econ');
        if ~isempty(D_max)
            [U,S,V,D_cut] = trim(U,S,V,D_max,epsilon);
            C = reshape(V',[D_cut,b_2,d]);
        else
            C = reshape(V',[size(V,2),b_2,d]);
        end
        carryover = U*S;
end
end

function [U,S,V,D_cut] = trim(U,S,V,D_max,epsilon)
% Computes appropriate trimming of U,S,V matrices after SVD
% The decimation is computed by taking the sum of the square of the singular
% values and then keeping only the ones under (1 - epsilon)
%
% INPUT
%   U,S,V:      matrices from SVD decompositon     
%   D_max:      maximum allowed bond length
%   epsilon:    tolerated error in decimation
% OUTPUT
%   U,S,V:      matrices from SVD decompositon after truncation
%   D_cut:      computed bond length

% Compute number of singular values under (1 - epsilon)
sum_s = cumsum(diag(S).^2);
D_cut = sum(sum_s < (1 - epsilon)*sum_s(end)) + 1;
if D_cut > D_max
    D_cut = D_max;
end
% Do trimming
U = U(:,1:D_cut);
S = S(1:D_cut,1:D_cut);
V = V(:,1:D_cut);
end
