function operator = expand_mpo(mpo)
% expands MPS into QM state vector

N = length(mpo);
d = size(mpo{1},4);

operator = zeros(d^N,d^N);
c = cell(1,N);  % The following code expands the MPS, may take a lot of mem
[c{:}] = ndgrid(1:d);
combs = fliplr(cell2mat(cellfun(@(v)v(:),c,'UniformOutput',false)));
for pos1 = 1:d^N
    for pos2 = 1:d^N
        prod = 1;
        for i = 1:N
            prod = prod*mpo{i}(:,:,combs(pos1,i),combs(pos2,i));
        end
        operator(pos1,pos2) = prod;
    end
end

end

