function result = multiply(O_1,O_2)
% Multiplies two MPO elements together with representation 
% (bond,bond,physical,physical).
% 
% INPUT
%	O_1,O_2:	MPO element, a rank-4 tensor with index convention
%				(bond,bond,physical,physical)
% OUTPUT
%	result:		corresponding MPO element after multiplication

result = ncon({O_2,O_1},{[-1,-3,-5,1],[-2,-4,1,-6]});
s = size(result);
result = reshape(result,[s(1)*s(2),s(3)*s(4),s(5),s(6)]);
end
