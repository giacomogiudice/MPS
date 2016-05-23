function melem = matrixelement(mps_1,mpo,mps_2)
% Returns the matrix element < MPS_1 | MPO | MPS_2 > that runs in
% polynomial time.
% The inputs are cell arrays corresponding to MPS decompositions.
% Notice that no check on correct sizes is done.

N = length(mps_1);
d = size(mps_1{1},3);
melem = 1;

for i = 1:N
    sum = 0;
    for sigma_1= 1:d
        for sigma_2 = 1:d
            sum = sum + mpo{i}(:,:,sigma_1,sigma_2)*(mps_1{i}(:,:,sigma_1)'*melem*mps_2{i}(:,:,sigma_2));
        end
    end
    melem = sum;
end
end