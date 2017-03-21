function [mps_out,iter,distance] = sweep_iter(mps_in,mpo,mps_out,varargin)
% Computes a DMRG-style sweep on an MPS, applying some operator in MPO form
% and then doing canonization and decimation using the iterative method.
%
% INPUT
%   mps_in:             cell array corresponding to input MPS
%   mpo:                cell array corresponding to MPO
%   mps_out:            cell array corresponding to guess of output MPS
%                       (WARNING: must be left canonized!)
%   iter_max:           (optional) maximum number of iterations
%   tolerance:          (optional) tolerance in the difference between the 
%                       guess and the full state
%   max_storage_size:   (optional) maximum number of bytes available to 
%                       store the full application of the MPO to  the MPS.
%                       If the actual size is smaller that this limit, the
%                       product MPO*MPS will be cached
% OUTPUT
%   mps_out:            resulting MPS after computation in right canonization
%   iter:               number of iterations in optimization
%   distance:           distance between the compressed MPS and the full MPS,
%                       useful to estimate the error in compressing

% Default values
iter_max = 10;
tolerance = 1e-4;
max_storage_size = 1e8;

if ~isempty(varargin)
    iter_max = varargin{1};
    tolerance = varargin{2}; 
    if length(varargin) > 2
        max_storage_size = varargin{3};
    end
end

if iter_max <= 0
    iter = 0;
    return
end

N = length(mps_in);

if isempty(mpo)
    mult = @(W,M,k) M;
    % Create dummy MPO
    mpo = cell(1,N);
else
    if predict_product_size(mpo,mps_in) <= max_storage_size
        mps_full = cellfun(@(W,M) apply(W,M),mpo,mps_in,'UniformOutput',false);
        mult = @(W,M,k) mps_full{k};
    else
        mult = @(W,M,k) apply(W,M);
    end
end

% Initialize block storage
blocks = cell(1,N+1);
blocks{1} = 1;
blocks{N+1} = 1;
for site = N:(-1):1
    target = mult(mpo{site},mps_in{site},site);
    blocks{site} = update_right(blocks{site+1},mps_out{site},target);
end

Kvalues = zeros(1,N);
for iter = 1:iter_max
    % Optimization sweep left -> right
    for site = 1:(N-1)
        % Apply MPO
        target = mult(mpo{site},mps_in{site},site);
        % Variational step
        mps_out{site} = optimization_step(target,blocks{site},blocks{site+1});
        % Canonize the new element
        mps_out{site} = canonize_fast(mps_out{site},+1);
        % Update current block
        blocks{site+1} = update_left(blocks{site},mps_out{site},target);
    end
    % Do same for last site, except update
    target = mult(mpo{N},mps_in{N},N);
    mps_out{N} = optimization_step(target,blocks{N},blocks{N+1});
    mps_out{N} = canonize_fast(mps_out{N},+1);
    % Optimization sweep right -> left
    for site = N:(-1):2
        % Apply MPO
        target = mult(mpo{site},mps_in{site},site);
        % Optimization step
        mps_out{site} = optimization_step(target,blocks{site},blocks{site+1});
        % calculate error
        Kvalues(site) = contract(conj(mps_out{site}),[1,2,3],mps_out{site},[1,2,3]);
        % Canonize the new element
        mps_out{site} = canonize_fast(mps_out{site},-1);
        % Update current block
        blocks{site} = update_right(blocks{site+1},mps_out{site},target);
    end
    % Do same for first site, except update
    target = mult(mpo{1},mps_in{1},1);
    mps_out{1} = optimization_step(target,blocks{1},blocks{2});
    Kvalues(1) = contract(conj(mps_out{1}),[1,2,3],mps_out{1},[1,2,3]);
    mps_out{1} = canonize_fast(mps_out{1},-1);
    % Overlap is the last update
    distance = 1 - abs(update_right(blocks{2},mps_out{1},target))
    % Calculate stopping condition
    if std(Kvalues)/abs(mean(Kvalues)) <= tolerance
        break;
    end
end
end

function new_guess = optimization_step(target,block_left,block_right)
    new_guess = contract(block_right,2,target,2);
    new_guess = contract(block_left,2,new_guess,2);
end

function new_block = update_left(prev_block,M_1,M_2)
    new_block = contract(prev_block,1,conj(M_1),1);
    new_block = contract(new_block,[1,3], M_2,[1,3]);
end

function new_block = update_right(prev_block,M_1,M_2)
    new_block = contract(prev_block,1,conj(M_1),2);
    new_block = contract(new_block,[1,3], M_2,[2,3]);
end

function bytes = predict_product_size(mpo,mps)
% Predicts the size of the resulting MPS in bytes if you compute the
% application of some MPO to an MPS provided. Notice that this is a 
% slightly optimistic estimate since it doesn't take in account the 
% MATLAB overhead.
% 
% INPUT
%   O:      1d cell array representing an MPO 
%   M:      1d cell array representing an MPS
% OUTPUT
%   nBytes: number of bytes corresponding to MPO*MPS

complexToBytes = 16;

nComplex = sum(cellfun(@(W,M) size(W,1)*size(W,2)*prod(size(M)),mpo,mps));
bytes = complexToBytes*nComplex;
end
