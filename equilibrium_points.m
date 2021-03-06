%% Initial Setup
clear; clc; close all;

% folders to create
image_folder = 'images/Equilibrium';
cache_folder = 'tmp/cache';

[~,~]=mkdir(image_folder);
image_folder = strcat(image_folder, '/');
[~,~]=mkdir(cache_folder);
cache_folder = strcat(cache_folder, '/');

addpath(genpath('functions'))
addpath(genpath('simulations'))
addpath(genpath('system_models'))
addpath(genpath('scripts'))

Simulink.fileGenControl('set', 'CacheFolder', cache_folder);

%% System Specifications

run system_specifications

%% Simulation Parametersr

% Desired DC-DC converter to use
% Options can be found in the system directory:
%   buck
%   boost
%   buck_boost
%   buck_boost_non_inverting
circuit = boost(R, Ro, Co, L);

%% Prepare Data

run load_circuit_sys

%% System Equilibrium Points

[sample_lambdas, equilibrium] = generate_sample_points(sys);

%% Analysis

plot_voltage_current_equilibrium(equilibrium, circuit.name, image_folder);
plot_voltage_lambda_equilibrium(equilibrium, sample_lambdas, circuit.name, image_folder);
