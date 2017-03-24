function X = contract(X,dimX,indX,Y,dimY,indY)
% Does tensor contraction between X and Y, summing over the indexes
% specified by indX and indY, in the order which they are stored in the 
% corresponding input tensors. MATLAB automatically discarting trailing
% singleton dimensions for tensors larger than rank-2. Because of this it
% is necessary to provide the dimension of the tensors X and Y. In this
% way one does not have to worry if the last indices are virtual and this
% will hopefully catch more coding mistakes.  
%
% INPUT
%   X, Y:       tensors to contract
%   indX, indY: array of indexes to contract of each tensor
%   dimX, dimY: an array of the number of indices in X and Y
% OUTPUT
%   X:          resulting tensor

% Compute the size of each tensor accounting for singleton dimensions
Xsize=ones(1,dimX); 
Xsize(1:ndims(X))=size(X);
Ysize=ones(1,dimY); 
Ysize(1:length(size(Y)))=size(Y);


% indXr and indYr contain the uncontracted indices
indXr=1:dimX; 
indXr(indX)=[];
indYr=1:dimY;
indYr(indY)=[]; 
% Define size of contracted and remaining indices
sizeXr=Xsize(indXr); 
sizeX=Xsize(indX); 
sizeYr=Ysize(indYr); 
sizeY=Ysize(indY); 
if any(sizeX ~= sizeY) % If any contracted dimension mismatch
    error('Dimension mismatch in indices provided');
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
