function mps_out = addMPS(mps_1,mps_2)
% Adds two MPS by combining each element in a block-diagonal way
% 
% INPUT
%   mps_1,mps_2:	MPS represented as a cell array, with 
%					index convention (bond,bond,physical,physical)
% OUTPUT
%   mps_2: 			MPS corresponding to the sum of the inputs

N = length(mps_1);
D = size(mps_1{1},3);
mps_out = cell(1,N);
% Do first element separately
for k = 1:D
	mps_out{1}(1,:,k) = [mps_1{1}(1,:,k),mps_2{1}(1,:,k)];
end
% Do the middle elements
for site = 2:(N-1)
	for k = 1:D
		mps_out{site}(:,:,k) = blkdiag(mps_1{site}(:,:,k),mps_2{site}(:,:,k));
	end
end
% Do last element separately
for k = 1:D
	mps_out{N}(:,1,k) = [mps_1{N}(:,1,k);mps_2{N}(:,1,k)];
end

end
