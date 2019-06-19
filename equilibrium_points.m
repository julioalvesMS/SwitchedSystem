%% Initial Setup
clear; clc; close all;

% folders to create
image_folder = 'images';
cache_folder = 'tmp/cache';

[~,~]=mkdir(image_folder);
image_folder = strcat(image_folder, '/');
[~,~]=mkdir(cache_folder);
cache_folder = strcat(cache_folder, '/');

addpath(genpath('functions'))
addpath(genpath('models'))

Simulink.fileGenControl('set', 'CacheFolder', cache_folder);

%% Simulation Parametersr

% Desired DC-DC converter to use
% Options can be found in the system directory:
%   buck
%   boost
%   buck_boost
circuit = buck;

sys = default_converter_sys(circuit);

%% System Equilibrium Points

[sample_lambdas, equilibrium] = generate_sample_points(sys);

%% Analysis

plot_voltage_current_equilibrium(equilibrium, circuit.name, image_folder);
plot_voltage_lambda_equilibrium(equilibrium, sample_lambdas, circuit.name, image_folder);
