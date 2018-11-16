function handle = proj_onesite(P,left,right)
M = length(P);
B = cell(M,1);
norms = zeros(M,1);
for m = 1:M
	B{m} = ncon({left{m},P{m},right{m}},{[-1,1],[1,2,-3],[-2,2]});
	norms(m) = B{m}(:)'*B{m}(:);
end
function W = new_element(M)
	W = M;
	for m = 1:m
		overlap = ncon({conj(B{m}),M},{[1,2,3],[1,2,3]});
		W = W - overlap/norms(m)*B{m};
	end
end
handle = @new_element;
end
