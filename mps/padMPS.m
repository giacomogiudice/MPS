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

bnd = zeros(2,N);
first_half = floor(N/2);
bnd(:,1:first_half) = [d.^((1:first_half)-1);d.^(1:first_half)];
if mod(N,2)
	first_half = first_half + 1;
	bnd(:,first_half) = [d^(first_half-1);d^(first_half-1)];	
end
first_half = first_half + 1;
bnd(:,first_half:N) = [d.^(N+1-(first_half:N));d.^(N-(first_half:N))];
bnd = min(bnd,D_max);

% Construct output MPS
for site = 1:N
	[D_1,D_2,~] = size(mps_in{site});
	if D_1 > D_max || D_2 > D_max
		disp('Unable to pad MPS, appears to be larger than size provided');
		mps_out = mps_in;
		return
	end
	diff_size = bnd(:,site)' - [D_1,D_2];
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
