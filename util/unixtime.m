function t = unixtime()
% Compute UNIX time, i.e. the number of seconds from January 1, 1970
%
% INPUT
%	(N/A)
% OUTPUT
%	t:	UNIX time as an integer
 
unix_origin = datenum('19700101 000000','yyyymmdd HHMMSS');
t = int32((now - unix_origin)*86400);
end
