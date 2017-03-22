% This is a simple script to test different MPS functions using the quantum
% Ising model in one dimension. It compares the MPS version with the 
% traditional matrix-vector operations. Expected output is that it passes 
% all assertions.
% Feel free to change N, but probably not the other parameters without
% knowing what you are doing.

close all
clear all
%% Parameters
N = 9;      % Number of sites
d = 2;      % Local H-space dimension
D_S = 2;    % State bond dimension
D_O = 3;    % Operator bond dimension

test_tolerance = 1e-12; % tolerance in numerical errors

%% MPS Initialization
M = cell(1,N);

A_odd = zeros(D_S,D_S,d);
A_odd(:,:,1) = eye(D_S);
A_even = zeros(D_S,D_S,d);
A_even(:,:,2) = eye(D_S);

M{1} = zeros(1,2,2);
M{1}(:,:,1) = [1,0];   
M{1}(:,:,2) = [0,0];

for i = 2:2:N-1
    M{i} = A_even;
    M{i + 1} = A_odd;
end

M{N} = zeros(2,1,2);
if mod(N,2) == 1
    M{N}(:,:,1) = [1;0];    
    M{N}(:,:,2) = [0;0];
else
    M{N}(:,:,1) = [0;0];    
    M{N}(:,:,2) = [1;0];
end

%% MPO Initialization
s_z = 1/sqrt(2)*[1,0;0,-1];
W = zeros(D_O,D_O,d,d);

W(1,1,:,:) = eye(2);
W(2,1,:,:) = s_z;
W(3,2,:,:) = s_z;
W(3,3,:,:) = eye(2);

O = cell(1,N);
O{1} = W(D_O,:,:,:);
for i = 2:N-1
    O{i} = W;
end
O{N} = W(:,1,:,:);

%% Traditional matrix Operations
neel = [1;0];
for i = 2:N
    if mod(i,2) == 0;
        neel = kron(neel,[0;1]);
    else
        neel = kron(neel,[1;0]);
    end
end

ham = zeros(2^N);
for i = 1:(N-1) % number of terms in H
    S = sparse(1);
    for j = 1:N
        if j == i
            S = kron(S,s_z);
        elseif j == i + 1
            S = kron(S,s_z);
        else
            S = kron(S,speye(2));
        end
    end
    ham = ham + S;
end

%% Calculate New State
Mprime = cell(1,N);
for i = 1:N
    Mprime{i} = apply(O{i},M{i});
end
newstate = expandMPS(Mprime);


%% Check Consistency with Vector Operations
fprintf('Testing MPS expansion...');
assert(norm(expandMPS(M) - neel) < test_tolerance);
fprintf('\tdone\n');
fprintf('Testing MPO expansion...');
assert(norm(expandMPO(O) - ham) < test_tolerance);
fprintf('\tdone\n');
fprintf('Testing MPO times MPS...');
assert(norm(newstate - ham*neel) < test_tolerance)
fprintf('\tdone\n');
fprintf('Testing scalar product...');
assert(norm(braket(M,Mprime) - neel'*ham*neel) < test_tolerance)
fprintf('\tdone\n');

%% Test Canonization with SVD
[Mleft,Mprime_norm] = sweep(M,O,+1);

fprintf('Testing left canonization...');
assert(iscanonized(Mleft,1,test_tolerance));
assert(1 - norm(braket(Mleft,Mprime)/Mprime_norm) < test_tolerance);
assert(1 - norm(expandMPS(Mleft)'*newstate/norm(newstate)) < test_tolerance);
fprintf('\tdone\n');

Mright = sweep(M,O,-1); 
fprintf('Testing right canonization...');
assert(iscanonized(Mright,-1,test_tolerance));
assert(1 - norm(braket(Mright,Mprime)/Mprime_norm) < test_tolerance);
assert(1 - norm(expandMPS(Mright)'*newstate/norm(newstate)) < test_tolerance);
fprintf('\tdone\n');

fprintf('All tests passed!\n');
