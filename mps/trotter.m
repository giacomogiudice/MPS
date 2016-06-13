function [U_odd,U_even] = trotter(U_pair,N,D)
% Creates MPO evolution operators for even and odd bonds to be used with
% Suzuki-Trotter decomposition
%
% INPUT
%   U_pair: evolution operator U = expm(-1i*H*dt) corresponding to the pair
%   N:      number of sites
%   D:      bond dimension
% OUTPUT
%   U_odd:  MPO evolution for odd bonds
%   U_even: MPO evolution for even bonds

%% U -> P reshaping and SVD, then more reshaping

P = reshape(U_pair,[D,D,D,D]);
P = permute(P,[1,3,2,4]);           % not P = permute(P,[4,2,3,1]);
P = reshape(P,[D^2,D^2]);           % (sigma_1 sigma_1'),(sigma_2 sigma_2')

[U,S,V] = svd(P,'econ');
k = size(S,2);
U = U*sqrt(S);                      % (sigma_1 sigma_1'), k
U_bar = sqrt(S)*V';                 % k,(sigma_2 sigma_2')

U = reshape(U,[D,D,1,k]);           % sigma_1,sigma_1',1,k
U = permute(U,[3,4,1,2]);           % 1,k,sigma_1,sigma_1'
U_bar = reshape(U_bar,[k,1,D,D]);   % k,1,sigma_2 sigma_2'
I = reshape(eye(D),[1,1,D,D]);      % 1,1,sigma,sigma'

%% Putting it all into U_even and U_odd
U_odd = cell(1,N);
U_even = cell(1,N);
for i = 1:2:N-1;
    U_odd{i}= U;
    U_odd{i+1}= U_bar;
    U_even{i} = U_bar;
    U_even{i+1} = U;
    
end
U_even{1}= I;

if mod(N,2)
    U_odd{N} = I; 
    U_even{N} = U_bar;
else
    U_even{N} = I;
end
end
