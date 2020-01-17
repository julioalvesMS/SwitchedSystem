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
addpath(genpath('simulations'))
addpath(genpath('system_models'))
addpath(genpath('scripts'))


plot_compression_rate = 1e0;


Simulink.fileGenControl('set', 'CacheFolder', cache_folder);

%% System Specifications

run system_specifications

%% Sistema

% Desired Theorem to use
% Theorems defines
%   1 - Fixed Equilibrium
%   2 - Valiable Equilibrium
opt_theorem = 1;

% Disturbances to be applied during simulations
% Options
%   disturbance_Vin_enable - Enable step disturbance in the input voltage
%   disturbance_Ro_enable - Enable step disturbance in the load resistance
disturbance_Vin_enable = false;
disturbance_Ro_enable = false;

circuit = boost(R, Ro, Co, L);

Vref = circuit.single_voltage;

run load_circuit_sys

simulation_duration = 1;

%% Generate Limit Cycle

xe = [0; Vref];
Gamma = [0 0.5];

[cand, kappa] = find_cycles(dsys, xe, Gamma);

opt = 1;
cycle = find_limit_cycle(dsys, kappa, cand, opt);

%% Prepare Data

run circuit_disturbance


%% Simulate Converter 

% Get model to simulate
model = circuit.limit_cycle_simulink;

load_system(model);

% Convert the truct used to represent the space state to double, so it
% can be used in the simulink
[Ad, Bd, Cd, Dd, Qd, Ld] = gss2double(dsys);    % Discrete

ell = zeros(length(sys.A{1}), kappa, sys.N);
for i=1:dsys.N
    ell(:,:,i) = cell2mat({cycle.ell{i,:}});
end

P = zeros(length(sys.A{1}), length(sys.A{1}), kappa);
for i=1:kappa
    P(:,:,i) = cycle.lyap.P{i};
end

xe_h = cell2mat(cycle.xe_h);

% Creates a bus, wich will be used in the simulink to simplify the
% model
SystemDataBus = create_bus_LimitCycleDataBus(Ad, Bd, P, Qd, sys.N, xe_h, ell);

% Run simulation
sim(model, simulation_duration);

i = 1;
% Store only samples of the data, this will be made in order to save
% memory use
sim_out(i).IL = downsample_timeseries(logsout.get('IL').Values, plot_compression_rate);
sim_out(i).Vout = downsample_timeseries(logsout.get('Vout').Values, plot_compression_rate);


%% Analysis

plot_voltage_current(sim_out, circuit.name, image_folder);

plot_voltage_time(sim_out, circuit.name, image_folder);

plot_current_time(sim_out, circuit.name, image_folder);

if disturbance_Ro_enable == 1
    plot_disturbance_voltage_time(sim_out, disturbance_Ro_time, circuit.name, image_folder);
end

if disturbance_Vin_enable == 1
    plot_disturbance_voltage_time(sim_out, disturbance_Vin_time, circuit.name, image_folder);
end