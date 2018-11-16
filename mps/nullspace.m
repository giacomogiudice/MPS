function [nulltensor,C,carryover] = nullspace(M,direction)
[D1,D2,d] = size(M);
if direction == 1
	Ddiff = D1*d - D2;
	A_left = reshape(permute(M,[1,3,2]),[d*D1,D2]);
	[Q,R] = qr(A_left);
	nulltensor = permute(reshape(Q(:,(D2+1):d*D1),[D1,d,Ddiff]),[1,3,2]);
	C = permute(reshape(Q(:,1:D2),[D1,d,D2]),[1,3,2]);
	carryover = R(1:D2,1:D2);
	F = ncon({C,carryover},{[-1,1,-3],[1,-2]}) - M;
	norm(F(:))
	G = ncon({conj(C),C},{[1,-1,2],[1,-2,2]}) - eye(D2);
	norm(G(:))
elseif direction == -1
	Ddiff = D2*d - D1;
	A_right = reshape(M,[D1,d*D2]);
	[Q,R] = qr(A_right.');
	nulltensor = reshape(Q(:,(D1+1):d*D2).',[Ddiff,D2,d]);
	C = reshape(Q(:,1:D1).',[D1,D2,d]);
	carryover = R(1:D1,1:D1).';
	F = ncon({carryover,C},{[-1,1],[1,-2,-3]) - M;
	norm(F(:))
	G = ncon({conj(C),C},{[-1,1,2],[-2,1,2]}) - eye(D2);
	norm(G(:))

else
	error(['Unrecognized direction ' direction]);
end
end
