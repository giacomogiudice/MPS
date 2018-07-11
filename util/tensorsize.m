function str = tensorsize(m)
% Returns a string with the size of the input tensor m. If the input is a
% cell array, it returns the size of each element in the cell array.
%
% INPUT	
%	m:		tensor or cell array of tensors
% OUTPUT
%	str:	string corresponding to the size of m, or each entry in m

sizetostr = @(s) strjoin(strtrim(cellstr(num2str(s.'))),'x');
if iscell(m)
	c = cellfun(@(x) sizetostr(size(x)),m,'UniformOutput',false);
	str = '';
	for ind = 1:size(c,1)
		str = sprintf('%s    [%s]\n',str,strjoin(c(ind,:),']    ['));
	end
else
	str = sprintf('    [%s]\n',sizetostr(size(m)));
end
end
