function Z = contract(X,dimX,indX,Y,dimY,indY)
% Does tensor contraction between X and Y, summing over the indexes
% specified by indX and indY, in the order which they are stored in the 
% corresponding input tensors. MATLAB automatically discarting trailing
% singleton dimensions for tensors larger than rank-2. Because of this it
% is necessary to provide the dimension of the tensors X and Y. In this
% way one does not have to worry if the last indices are virtual and this
% will hopefully catch more coding mistakes.  
% The contraction is achieved by reshuffling the indices of the tensors 
% so that all the indices to contract are placed at the end of X and at 
% the beginning of Y. In this way the problem can be reduced to a matrix
% multiplication, which MATLAB is very fast at.
%
% INPUT
%	X, Y:	tensors to contract
%	indX, indY:	array of indexes to contract of each tensor
%	dimX, dimY:	an array of the number of indices in X and Y
% OUTPUT
%	Z:			resulting tensor after contracting (X,Y) over (indX,indY)

% Compute the size of each tensor accounting for singleton dimensions
Xsize=ones(1,dimX); 
Xsize(1:ndims(X))=size(X);
Ysize=ones(1,dimY); 
Ysize(1:length(size(Y)))=size(Y);

% indXr and indYr contain the remaining (uncontracted) indices
indXr=1:dimX; 
indXr(indX)=[];
indYr=1:dimY;
indYr(indY)=[]; 
% Define size of summed indices and remaining indices
sizeXr=Xsize(indXr); 
sizeX=Xsize(indX); 
sizeYr=Ysize(indYr); 
sizeY=Ysize(indY); 
% If any contracted dimension mismatch
if any(sizeX ~= sizeY)
	error('Dimension mismatch in indices provided');
end

if isempty(indYr)
	% If both tensors are fully contracted, place the indices in the correct
	% order, reshape both into a vector and perform the contraction as a 
	% vector-vector multiplication
	if isempty(indXr) 
		X=permute(X,indX);
		X=reshape(X,[1,prod(sizeX)]); 
		Y=permute(Y,indY); 
		Y=reshape(Y,[prod(sizeY),1]);
		Z=X*Y; 
		return
	% If Y is fully contracted, place the indices in the correct order,
	% reshape Y into a vector and perform a matrix-vector operation
	else 
		X=permute(X,[indXr,indX]);
		X=reshape(X,[prod(sizeXr),prod(sizeX)]);
		Y=permute(Y,indY);
		Y=reshape(Y,[prod(sizeY),1]);
		Z=X*Y;
		Xsize=Xsize(indXr);
		Z=reshape(Z,[Xsize,1]);
		return
	end
end
% Otherwise, both have uncontracted indices. Send summed indices to the 
% last places of X, send the summed indices to the first places of Y, and
% perform the contraction as a matrix-matrix operation
X=permute(X,[indXr,indX]);
X=reshape(X,[prod(sizeXr),prod(sizeX)]);
Y=permute(Y,[indY,indYr]);
Y=reshape(Y,[prod(sizeY),prod(sizeYr)]);
Z=X*Y;
Z=reshape(Z,[Xsize(indXr),Ysize(indYr)]);
end
