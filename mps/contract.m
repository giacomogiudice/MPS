function X = contract(X,indX,Y,indY)
% Does tensor contraction between X and Y, summing over the indexes
% specified by indX and inY
%
% INPUT
%   X,Y: tensors to contract
%   indX,indY: array of indexes to contract of each tensor
% OUTPUT
%   X: resulting tensor

numindX = length(size(X)); % total number of indexes for each tensor
numindY = length(size(Y)); 

Xsize=size(X); % Size of input tensors
Ysize=size(Y);

indXr=1:numindX; 
indXr(indX)=[];
indYr=1:numindY;
indYr(indY)=[]; % indXr and indYr are vectors of uncontracted indices
sizeXr=Xsize(indXr); 
sizeX=Xsize(indX); 
sizeYr=Ysize(indYr); 
sizeY=Ysize(indY); % size of contracted and remaining indices
if any(sizeX ~= sizeY) % If any contracted dimension mismatch
    error('indX and indY are not of same dimension.');
end
if isempty(indYr)
    if isempty(indXr) % If no uncontracted indexes for both
        X=permute(X,indX);
        X=reshape(X,[1,prod(sizeX)]);
        % Arrange them by correct contraction indexes, then reshape into
        %vector
        Y=permute(Y,indY);
        Y=reshape(Y,[prod(sizeY),1]); 
        X=X*Y; % Output matrix
        return
    else % If no uncontracted indexes for Y but not for X
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

% If both have uncontracted indices
X=permute(X,[indXr,indX]); % Concatenated array indXl,indX. Sends summed indexes to last places
X=reshape(X,[prod(sizeXr),prod(sizeX)]); % Matrix reshape
Y=permute(Y,[indY,indYr]); % Sends summed to first places
Y=reshape(Y,[prod(sizeY),prod(sizeYr)]); % Matrix reshape
X=X*Y;
Xsize=[Xsize(indXr),Ysize(indYr)]; % Size vector of final tensor
X=reshape(X,Xsize); % Reshape to final size