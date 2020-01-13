function H = hamiltonian(N,varargin)
% Build a Hamiltonian from local interactions. Each input is expected to
% be a cell array of the form {O_1,O_2,...} or {c,O_1,...}, where c is a
% scalar and O_n is an operator acting on a single site.

flag_pbc = false;
flag_reduce = true;
pbc = @(n) mod(n+N-1,N) + 1;

% n_terms = length(varargin);
j = 1;
% Preprocessing step
while j <= length(varargin)
    t = varargin{j};
    if ischar(t)
        if isequal(t,'-pbc');
            flag_pbc = true;
        elseif isequal(t,'-noreduce')
            flag_reduce = false;
        end
        varargin(j) = [];
    elseif iscell(t)
        if isscalar(t{1})
            % Throw the scalar in the first operator
            coeff = t{1};
            varargin{j}(1) = [];
            varargin{j}{1} = coeff*varargin{j}{1};
        else
            d = size(t{1},1);
            assert(all(cellfun(@(t) isequal(size(t),[d,d]),t)),'All operators are expected to be same size');
        end
        j = j + 1;
    end
end

chi = sum(cellfun(@(t) length(t) - 1,varargin)) + 2;
if flag_pbc
    n_pbc_terms = sum(cellfun(@(t) length(t) > 1,varargin));
else
    n_pbc_terms = 0;
end
W = zeros([chi+n_pbc_terms,chi+n_pbc_terms,d,d]);
I = reshape(eye(d),[1,1,d,d]);
W(1,1,:,:) = I;
W(chi,chi,:,:) = I;
if flag_pbc
    for k = 1:n_pbc_terms
        W(chi+k,chi+k,:,:) = I;
    end
end

% Now populate W in the upper-triangular part
ind = 1;
for j = 1:length(varargin)
    t = varargin{j};
    if length(t) == 1
        % Single-body term
        O = reshape(t{1},[1,1,d,d]);
        W(1,chi,:,:) = W(1,chi,:,:) + O;
    else
        for k = 1:length(t)
            O = reshape(t{k},[1,1,d,d]);
            if k == 1
                W(1,ind+1,:,:) = O;
            elseif k == length(t)
                W(ind+k-1,chi,:,:) = O;
            else
                W(ind+k-1,ind+k,:,:) = O;
            end
        end
        ind = ind + (k-1);
    end
end
% Build full operator
H = cell(1,N);
for site = 2:(N-1)
    H{site} = W;
end
H{1} = W(1,:,:,:);
H{N} = W(:,chi,:,:);

% Add PBC terms
if flag_pbc
    ind = 1;
    for j = 1:length(varargin)
        t = varargin{j};
        if length(t) > 1
            for k = 1:length(t)
                O = reshape(t{k},[1,1,d,d]);
                site = pbc(N+k-1);
                if site == 1
                    H{site}(1,chi+ind,:,:) = O;
                elseif site == N
                    H{site}(chi+ind,1,:,:) = O;
                else
                    H{site}(chi+ind,chi+ind,:,:) = O;
                end
            end
            ind = ind + 1;
        end
    end
end
if flag_reduce
    H = reduceMPO({H});
end
end
