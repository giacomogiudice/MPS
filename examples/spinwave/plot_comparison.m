% This script plots the comparison of the different methods for the same time
% step. Please run some of the simulations 'heisenberg_*' with the same
% parameters before comparing. The data of these simulations is expected
% to be stored in 'filename'.


close all;
setup;
load(filename);

%% Retrieve data
profiles = {};
legend_str = {};
if exist('magn_conventional','var')
	disp('Data for conventional evolution found...')
	profiles{end+1} = magn_conventional;
	legend_str{end+1} = 'Conventional';
end
if exist('magn_svd','var')
	disp('Data for SVD evolution found...')
	profiles{end+1} = magn_svd;
	legend_str{end+1} = 'ST2 (SVD)';
end
if exist('magn_iter','var')
	disp('Data for iterative evolution found...')
	profiles{end+1} = magn_iter;
	legend_str{end+1} = 'ST2 (iterative)';
end
if exist('magn_tdvp','var')
	disp('Data for TDVP found...')
	profiles{end+1} = magn_tdvp;
	legend_str{end+1} = 'TDVP';
end
if exist('magn_tdvp2','var')
	disp('Data for TDVP2 evolution found...')
	profiles{end+1} = magn_tdvp2;
	legend_str{end+1} = 'TDVP2';
end

if isempty(profiles)
	disp('Please run some simulations first!')
	return
end

figure
sites = 1:N;
hold on
for k = 1:length(profiles)
	handles(k) = plot(sites,profiles{k}(:,1));
end
legend(legend_str);
xlabel('sites')
ylabel('magnetization')
axis([0.5 N+0.5 -1.2 1.2])
for step = 1:time_steps
	pause(1/16)
	for k = 1:length(profiles)
		set(handles(k),'YData',profiles{k}(:,step));
	end
end
