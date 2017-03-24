function result = apply(O,M)
% Applies an MPO element O to an MPS element M
% The MPS index representation is assumed to be (bond,bond,physical) and
% the corresponding MPO representation should be 
% (bond,bond,physical,physical).
% 
% INPUT
%   O:      MPO element, a rank-4 tensor with index convention
%           (bond,bond,physical,physical)
%   M:      MPS element, a rank-3 tensor with index convention 
%           (bond,bond,physical)
% OUTPUT
%   result: corresponding MPS element after application of MPO

result = contract(M,3,3,O,4,4);
s = size(result);
result = reshape(permute(result,[3,1,4,2,5]),[s(1)*s(3),s(2)*s(4),s(5)]);
end
