function [C,carryover] = canonize(M,direction,varargin)
% Computes canonization of MPS site M according
% to left- (direction == +1) or right- canonization (direction == -1).
% It then optionally performs a decimation by a certain bond size if
% specified by the last two arguments
%
% INPUT
%   M:          site element of MPS, represented as a 3D tensor
%   direction:  specifies left (-1) or right (+1) canonization
%   D_max:      (optional) maximum allowed bond length
%   epsilon:    (optional) tolerated error in decimation
% OUTPUT
%   C:          MPS element canonized, A or B
%   carryover:  product to multiply with next neighbor, S*V' or U*S

if isempty(varargin)
    D_max = [];
    epsilon = [];
else
    D_max = varargin{1};
    epsilon = varargin{2};
end
s_M = size(M);

switch direction
    case +1 % Going right
        C = reshape(permute(M, [1,3,2]),s_M(3)*s_M(1),s_M(2));
        [U,S,V] = svd(C,'econ');
        if ~isempty(D_max)
            trim;
            % reshape V'
            C = permute(reshape(U,[s_M(1),s_M(3),D_cut]),[1 3 2]);
        else
            C = permute(reshape(U,[s_M(1),s_M(3),size(U,2)]),[1 3 2]);
        end
        carryover = S*V';
        
    case -1 % Going left
        C = reshape(M,[s_M(1),s_M(2)*s_M(3)]);
        [U,S,V] = svd(C,'econ');
        if ~isempty(D_max)
            trim;
            % reshape V'
            C = reshape(V',[D_cut,s_M(2),s_M(3)]);
        else
            C = reshape(V',[size(V,2),s_M(2),s_M(3)]);
        end
        carryover = U*S;
end
    function trim
        % Compute number of singular values under (1 - epsilon)
        sum_s = cumsum(gather(diag(S)).^2);
        D_cut = sum(sum_s < (1 - epsilon)*sum_s(end)) + 1;
        if D_cut > D_max
            D_cut = D_max;
        end
        % Do trimming
        U = U(:,1:D_cut);
        S = S(1:D_cut,1:D_cut);
        S = S/norm(diag(S));
        V = V(:,1:D_cut);
    end
end
