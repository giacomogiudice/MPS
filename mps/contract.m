function X = contract(X,indX,Y,indY)
% Does tensor contraction between X and Y, summing over the indexes
% specified by indX and inY, in the order which they are stored in the 
% corresponding input tensors.
%
% INPUT
%   X, Y:       tensors to contract
%   indX, indY: array of indexes to contract of each tensor
% OUTPUT
%   X:          resulting tensor

% Total number of indices for each tensor
numindX = length(size(X)); 
numindY = length(size(Y)); 
% Size of input tensors
Xsize=size(X); 
Ysize=size(Y);
% indXr and indYr are vectors of uncontracted indices
indXr=1:numindX; 
indXr(indX)=[];
indYr=1:numindY;
indYr(indY)=[]; 
% Define size of contracted and remaining indices
sizeXr=Xsize(indXr); 
sizeX=Xsize(indX); 
sizeYr=Ysize(indYr); 
sizeY=Ysize(indY); 
if any(sizeX ~= sizeY) % If any contracted dimension mismatch
    error('indX and indY are not of same dimension.');
end
if isempty(indYr)
    if isempty(indXr) % If no uncontracted indices for both
        X=permute(X,indX);
        X=reshape(X,[1,prod(sizeX)]); 
        Y=permute(Y,indY); % Arrange X and Y by correct contraction indexes
        Y=reshape(Y,[prod(sizeY),1]); % Reshape into vector
        X=X*Y; % Output matrix
        return
    else % If no uncontracted indices for Y but not for X
        X=permute(X,[indXr,indX]);
        X=reshape(X,[prod(sizeXr),prod(sizeX)]);
        Y=permute(Y,indY);
        Y=reshape(Y,[prod(sizeY),1]);
        X=X*Y;
        Xsize=Xsize(indXr);
        X=reshape(X,[Xsize,1]);
        return
    end
end
% Otherwise,  both have uncontracted indices
X=permute(X,[indXr,indX]); %  Send summed indices to last places
X=reshape(X,[prod(sizeXr),prod(sizeX)]); % Matrix reshape
Y=permute(Y,[indY,indYr]); % Send summed indices to first places
Y=reshape(Y,[prod(sizeY),prod(sizeYr)]); % Matrix reshape
X=X*Y; % Perform contraction as a matrix-matrix operation
X=reshape(X,[Xsize(indXr),Ysize(indYr)]); % Reshape to final size
end
