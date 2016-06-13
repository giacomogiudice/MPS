function result = apply(O,M)
% Applies an MPO element O to an MPS element M
% The MPS index representation is assumed to be (bond,bond,physical) and
% the corresponding MPO representation should be
% 
% INPUT
%   O:      MPO element, a rank-4 tensor with index convention
%           (bond,bond,physical,physical)
%   M:      MPS element, a rank-3 tensor with index convention 
%           (bond,bond,physical)
% OUTPUT
%   result: corresponding MPS element after application of MPO

result = contract(M,3,O,4);
s_r = size(result);
result = reshape(permute(result,[3 1 4 2 5]),[s_r(1)*s_r(3),s_r(2)*s_r(4),s_r(5)]);

end