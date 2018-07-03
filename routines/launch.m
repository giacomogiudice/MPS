function launch(filename,varargin)
% Used to launch a simulation from the terminal. Save the final state
% to a file named after the input filename
%
% INPUT
%   filename:   .mat file of input parameters
%   varargin:   pair(s) of variable names and values that override those
%               loaded from the file

if mod(length(varargin),2) ~= 0
    fprintf('Error: expecting even number of parameters, have %d\n',length(varargin));
    return
end

load(filename);

showbar = 0;
fileout = strcat(filename,'__');
for k = 1:2:length(varargin)
    eval([varargin{k} '=' varargin{k+1} ';']);
    fprintf('Setting %s = %s\n',varargin{k},varargin{k+1});
    fileout = strcat(fileout,varargin{k},varargin{k+1});
end
clear filename varargin;

[~,parrank] = system('echo $SLURM_PROCID');
parrank = str2num(parrank);
if ~isempty(parrank) 
    rng(unixtime + parrank);
    fileout = [fileout '_r' num2str(parrank)];
else
    parrank = 0;
    disp('Warning: no parallel rank detected');
end
fprintf('Saving at "%s". Running...\n',fileout)

tic
disp('<Insert your routine here>')
%    ...    %
toc
fprintf('Rank %d done\n',parrank);
end
