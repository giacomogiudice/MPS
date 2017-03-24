function result = multiply(O_1,O_2)
% Multiplies two MPO elements together with representation 
% (bond,bond,physical,physical).
% 
% INPUT
%   O_1,O_2:    MPO element, a rank-4 tensor with index convention
%				(bond,bond,physical,physical)
% OUTPUT
%   result:     corresponding MPO element after multiplication

result = contract(O_2,4,3,O_1,4,4);
s = size(result);
result = reshape(permute(result,[4,1,5,2,6,3]),[s(1)*s(4),s(2)*s(5),s(6),s(3)]);
end
