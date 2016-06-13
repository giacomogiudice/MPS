function launch_simulation(filename,varargin)
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
fileout = strcat(filename,'_out');
for i = 1:2:length(varargin)
    eval([varargin{i} '=' varargin{i+1} ';']);
    fprintf('Setting %s = %s\n',varargin{i},varargin{i+1});
    fileout = strcat(fileout,varargin{i},varargin{i+1});
end
fileout = strcat(fileout,'.mat');
fprintf('Running...\n')
tic
density_mps;
toc
fprintf('Complete. Saving as "%s"\n',fileout)
save(fileout)
fprintf('Done\n')
end
