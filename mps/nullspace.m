function [nulltensor,C,carryover] = nullspace(M,direction)
[D1,D2,d] = size(M);
if direction == 1
	Ddiff = D1*d - D2;
	A_left = reshape(permute(M,[1,3,2]),[d*D1,D2]);
	[Q,R] = qr(A_left);
	nulltensor = permute(reshape(Q(:,(D2+1):d*D1),[D1,d,Ddiff]),[1,3,2]);
	C = permute(reshape(Q(:,1:D2),[D1,d,D2]),[1,3,2]);
	carryover = R(1:D2,1:D2);
elseif direction == -1
	Ddiff = D2*d - D1;
	A_right = reshape(M,[D1,d*D2]);
	[Q,R] = qr(A_right.');
	nulltensor = reshape(Q(:,(D1+1):d*D2).',[Ddiff,D2,d]);
	C = reshape(Q(:,1:D1).',[D1,D2,d]);
	carryover = R(1:D1,1:D1).';
else
	error(['Unrecognized direction ' direction]);
end
end
