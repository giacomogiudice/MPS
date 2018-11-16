function result = apply(O,M)
% Applies an MPO element O to an MPS element M
% The MPS index representation is assumed to be (bond,bond,physical) and
% the corresponding MPO representation should be 
% (bond,bond,physical,physical).
% 
% INPUT
%	O:	MPO element, a rank-4 tensor with index convention
%		(bond,bond,physical,physical)
%	M:	MPS element, a rank-3 tensor with index convention 
%		(bond,bond,physical)
% OUTPUT
%	result:	corresponding MPS element after application of MPO

result = ncon({O,M},{[-1,-3,-5,1],[-2,-4,1]});
s = size(result);
result = reshape(result,[s(1)*s(2),s(3)*s(4),s(5)]);
end
