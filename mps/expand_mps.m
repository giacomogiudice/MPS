function state = expand_mps(mps)
% expands MPS into QM state vector

N = length(mps);
d = size(mps{1},3);

state = zeros(d^N,1);
c = cell(1,N);
[c{:}] = ndgrid(1:d);
combs = fliplr(cell2mat(cellfun(@(v)v(:),c,'UniformOutput',false)));
for pos = 1:d^N
    prod = 1;
    for i = 1:N
        prod = prod*mps{i}(:,:,combs(pos,i));
    end
    state(pos) = prod;
end

end

