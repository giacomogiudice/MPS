% This script initializes all the parameters necessary for the simulation
% of driven-dissipative Bose-Hubbard chains. After this script is called,
% any simulation script can be called.

%% Wipe Everything
clear 
close all
clc

addpath('util','mps','routines','workspaces');

%% Default Simulation Parameters
N = 7;              % Number of sites
d = 3;              % Number of boson energy levels
D_max = 100;        % Maximum bond dimension
epsilon = 1e-12;    % Error tolerated in truncation
U = 5;              % Interaction energy
J = 2;              % Hopping energy
F = 1;              % Driving energy
d_omega = -2*J;     % Frequency difference omega_0 - omega_p
gamma = 1;          % Dissipation rate

T = 20;             % Stopping time
dt = 1e-2;          % Time step
n_sampling = 1e1;   % Sampling rate for measurements (measurements every 
                    % n_sampling*dt)

state = 0;          % Is interpreted as vacuum by the main routines
showbar = 1;        % If true shows interactive waitbar

%% Save to Workspace
save('workspaces/default.mat')
