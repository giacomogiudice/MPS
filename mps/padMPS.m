function mps_out = padMPS(mps_in,D_max,varargin)
% Pads an MPS with zeros, increasing bond indices to the largest possible
% value for open boundary conditions (OBC) up to D_max specified. Optionally 
% canonized the outputted state by specifying a direction.
%
% INPUT
%	mps_in:		cell array corresponding to input MPS
%	D_max:		maximum bond of output MPS
%	direction:	(optional) direction of canonization for output MPS
% OUTPUT
%	mps_out:	resulting MPS with additional entries of zeros

N = length(mps_in);
d = size(mps_in{1},3);
mps_out = cell(1,N);
full_size = cell(1,N);

% Build cell array containing size of output MPS
for site = 1:(N/2)
	full_size{site} = [min(d^(site-1),D_max),min(d^site,D_max)];
end
if mod(N,2)
	site = site + 1;
	full_size{site} = [min(d^(site-1),D_max),min(d^(site-1),D_max)];
end
for site = (site+1):N
	full_size{site} = [min(d^(N-site+1),D_max),min(d^(N-site),D_max)];
end

% Construct output MPS
for site = 1:N
	[D_1,D_2,~] = size(mps_in{site});
	if D_1 > D_max || D_2 > D_max
		disp('Unable to pad MPS, appears to be larger than size provided');
		mps_out = mps_in;
		return
	end
	diff_size = full_size{site} - [D_1,D_2];
	if any(diff_size > 0)
		for k = 1:d
			mps_out{site}(:,:,k) = blkdiag(mps_in{site}(:,:,k),zeros(diff_size));
		end
	else
		mps_out{site} = mps_in{site};
	end
end

% If necessary, canonize output MPS
if ~isempty(varargin)
	if varargin{1} == 1 || varargin{1} == -1
		mps_out = sweep_fast(mps_out,varargin{1});
	end
end
end
