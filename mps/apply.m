function result = apply(O,M)
% Applies an MPO element O to an MPS element M

result = contract(M,3,O,4);
s_r = size(result);
result = reshape(permute(result,[3 1 4 2 5]),[s_r(1)*s_r(3),s_r(2)*s_r(4),s_r(5)]);

end