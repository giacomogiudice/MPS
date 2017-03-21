function [U,U_bar] = trotter(U_pair,D)
% Creates MPO evolution elements for a pair of bonds to be used with
% Suzuki-Trotter decomposition.
%
% INPUT
%   U_pair: evolution operator U = expm(-1i*H*dt) corresponding to the pair
%   D:      bond dimension
% OUTPUT
%   U:  	MPO evolution element for first site
%   U_bar: 	MPO evolution element for second

% U -> P reshaping and SVD, then more reshaping
P = reshape(U_pair,[D,D,D,D]);
P = permute(P,[4,2,3,1]);
% P = permute(P,[1,3,2,4]);           % not P = permute(P,[4,2,3,1]);
P = reshape(P,[D^2,D^2]);           % (sigma_1 sigma_1'),(sigma_2 sigma_2')

[U,S,V] = svd(P,'econ');
k = size(S,2);
U = U*sqrt(S);                      % (sigma_1 sigma_1'), k
U_bar = sqrt(S)*V';                 % k,(sigma_2 sigma_2')

U = reshape(U,[D,D,1,k]);           % sigma_1,sigma_1',1,k
U = permute(U,[3,4,1,2]);           % 1,k,sigma_1,sigma_1'
U_bar = reshape(U_bar,[k,1,D,D]);   % k,1,sigma_2 sigma_2'
end
