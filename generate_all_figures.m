%% Initial Setup
clear; clc; close all;

% folders to create
root_image_folder = 'images';
cache_folder = 'tmp/cache';

[~,~]=mkdir(root_image_folder);
root_image_folder = strcat(root_image_folder, '/');
[~,~]=mkdir(cache_folder);
cache_folder = strcat(cache_folder, '/');

addpath(genpath('functions'))
addpath(genpath('simulations'))
addpath(genpath('system_models'))
addpath(genpath('scripts'))


plot_compression_rate = 1e3;


Simulink.fileGenControl('set', 'CacheFolder', cache_folder);

%% System Specifications

run system_specifications

run circuit_disturbance

%% Default Parameters

% Desired DC-DC converter to use
% Options can be found in the system directory:
%   buck
%   boost
%   buck_boost
circuit_buck = buck(R, Ro, Co, L);
circuit_boost = boost(R, Ro, Co, L);
circuit_buck_boost = buck_boost(R, Ro, Co, L);

config_simulations = {};

default_config.opt_model = 2;
default_config.opt_discrete = 0;
default_config.opt_theorem = 1;
default_config.opt_pwm = 0;
default_config.opt_equilibrium_controller = 0;
default_config.opt_update_equilibrium = 1;
default_config.opt_partial_information = 0;
default_config.opt_constant_reference = 1;
default_config.disturbance_Vin_enable = 0;
default_config.disturbance_Ro_enable = 0;
default_config.disturbance_Vin_time = disturbance_Vin_time;
default_config.disturbance_Ro_time = disturbance_Ro_time;
default_config.switching_period = -1;
default_config.circuit = [];
default_config.test_voltages = [];
default_config.simulation_duration = 0.15;
default_config.image_folder = root_image_folder;


%% Run PWM

image_folder = strcat(root_image_folder, '/PWM');
[~,~]=mkdir(image_folder);
image_folder = strcat(image_folder, '/');

pwm_config = default_config;
pwm_config.opt_pwm = 1;
pwm_config.image_folder = image_folder;
pwm_config.simulation_duration = 0.1;

% Buck
new_sim = pwm_config;
new_sim.circuit = circuit_buck;
new_sim.test_voltages = new_sim.circuit.test_voltages;
config_simulations{end+1} = new_sim;

% Boost
new_sim = pwm_config;
new_sim.circuit = circuit_boost;
new_sim.test_voltages = new_sim.circuit.test_voltages;
config_simulations{end+1} = new_sim;

% Buck-Boost
new_sim = pwm_config;
new_sim.circuit = circuit_buck_boost;
new_sim.test_voltages = new_sim.circuit.test_voltages;
config_simulations{end+1} = new_sim;

%% Run Theorem 1

image_folder = strcat(root_image_folder, '/Theorem_1');
[~,~]=mkdir(image_folder);
image_folder = strcat(image_folder, '/');

theo1_config = default_config;
theo1_config.opt_theorem = 1;
theo1_config.image_folder = image_folder;
theo1_config.simulation_duration = 0.1;

% Buck
new_sim = theo1_config;
new_sim.circuit = circuit_buck;
new_sim.test_voltages = new_sim.circuit.test_voltages;
config_simulations{end+1} = new_sim;

% Boost
new_sim = theo1_config;
new_sim.circuit = circuit_boost;
new_sim.test_voltages = new_sim.circuit.test_voltages;
config_simulations{end+1} = new_sim;

% Buck-Boost
new_sim = theo1_config;
new_sim.circuit = circuit_buck_boost;
new_sim.test_voltages = new_sim.circuit.test_voltages;
config_simulations{end+1} = new_sim;


%% Run Theorem 2

image_folder = strcat(root_image_folder, '/Theorem_2');
[~,~]=mkdir(image_folder);
image_folder = strcat(image_folder, '/');

theo2_config = default_config;
theo2_config.opt_theorem = 2;
theo2_config.image_folder = image_folder;

% Buck
new_sim = theo2_config;
new_sim.circuit = circuit_buck;
new_sim.test_voltages = new_sim.circuit.test_voltages;
new_sim.simulation_duration = 0.15;
config_simulations{end+1} = new_sim;

% Boost
new_sim = theo2_config;
new_sim.circuit = circuit_boost;
new_sim.test_voltages = new_sim.circuit.test_voltages;
new_sim.simulation_duration = 1.5;
config_simulations{end+1} = new_sim;

% Buck-Boost
new_sim = theo2_config;
new_sim.circuit = circuit_buck_boost;
new_sim.test_voltages = new_sim.circuit.test_voltages;
new_sim.simulation_duration = 1.5;
config_simulations{end+1} = new_sim;


%% Run Theorem 2 PI

image_folder = strcat(root_image_folder, '/Theorem_2_PI');
[~,~]=mkdir(image_folder);
image_folder = strcat(image_folder, '/');

theo2pi_config = default_config;
theo2pi_config.opt_theorem = 2;
theo2pi_config.opt_equilibrium_controller = 1;
theo2pi_config.image_folder = image_folder;

% Buck
new_sim = theo2pi_config;
new_sim.circuit = circuit_buck;
new_sim.test_voltages = new_sim.circuit.test_voltages;
new_sim.simulation_duration = 0.15;
config_simulations{end+1} = new_sim;

% Boost
new_sim = theo2pi_config;
new_sim.circuit = circuit_boost;
new_sim.test_voltages = new_sim.circuit.test_voltages;
new_sim.simulation_duration = 0.7;
config_simulations{end+1} = new_sim;

% Buck-Boost
new_sim = theo2pi_config;
new_sim.circuit = circuit_buck_boost;
new_sim.test_voltages = new_sim.circuit.test_voltages;
new_sim.simulation_duration = 0.7;
config_simulations{end+1} = new_sim;


%% Run PWM - Load Change

image_folder = strcat(root_image_folder, '/PWM/Load_Change');
[~,~]=mkdir(image_folder);
image_folder = strcat(image_folder, '/');

pwm_config = default_config;
pwm_config.opt_pwm = 1;
pwm_config.image_folder = image_folder;
pwm_config.disturbance_Ro_enable = 1;
pwm_config.simulation_duration = 0.6;

% Buck
new_sim = pwm_config;
new_sim.circuit = circuit_buck;
new_sim.test_voltages = new_sim.circuit.single_voltage;
config_simulations{end+1} = new_sim;

% Boost
new_sim = pwm_config;
new_sim.circuit = circuit_boost;
new_sim.test_voltages = new_sim.circuit.single_voltage;
config_simulations{end+1} = new_sim;

% Buck-Boost
new_sim = pwm_config;
new_sim.circuit = circuit_buck_boost;
new_sim.test_voltages = new_sim.circuit.single_voltage;
config_simulations{end+1} = new_sim;

%% Run Theorem 1 - Load Change

image_folder = strcat(root_image_folder, '/Theorem_1/Load_Change');
[~,~]=mkdir(image_folder);
image_folder = strcat(image_folder, '/');

theo1_config = default_config;
theo1_config.opt_theorem = 1;
theo1_config.image_folder = image_folder;
theo1_config.disturbance_Ro_enable = 1;
theo1_config.disturbance_Ro_time = 0.25;
theo1_config.simulation_duration = 0.45;

% Buck
new_sim = theo1_config;
new_sim.circuit = circuit_buck;
new_sim.test_voltages = new_sim.circuit.single_voltage;
config_simulations{end+1} = new_sim;

% Boost
new_sim = theo1_config;
new_sim.circuit = circuit_boost;
new_sim.test_voltages = new_sim.circuit.single_voltage;
config_simulations{end+1} = new_sim;

% Buck-Boost
new_sim = theo1_config;
new_sim.circuit = circuit_buck_boost;
new_sim.test_voltages = new_sim.circuit.single_voltage;
config_simulations{end+1} = new_sim;


%% Run Theorem 2 - Load Change

image_folder = strcat(root_image_folder, '/Theorem_2/Load_Change');
[~,~]=mkdir(image_folder);
image_folder = strcat(image_folder, '/');

theo2_config = default_config;
theo2_config.opt_theorem = 2;
theo2_config.image_folder = image_folder;
theo2_config.disturbance_Ro_enable = 1;

% Buck
new_sim = theo2_config;
new_sim.circuit = circuit_buck;
new_sim.test_voltages = new_sim.circuit.single_voltage;
new_sim.disturbance_Ro_time = 0.2;
new_sim.simulation_duration = 0.45;
config_simulations{end+1} = new_sim;

% Boost
new_sim = theo2_config;
new_sim.circuit = circuit_boost;
new_sim.test_voltages = new_sim.circuit.single_voltage;
new_sim.disturbance_Ro_time = 1.8;
new_sim.simulation_duration = 2;
config_simulations{end+1} = new_sim;

% Buck-Boost
new_sim = theo2_config;
new_sim.circuit = circuit_buck_boost;
new_sim.test_voltages = new_sim.circuit.single_voltage;
new_sim.disturbance_Ro_time = 1.8;
new_sim.simulation_duration = 2;
config_simulations{end+1} = new_sim;


%% Run Theorem 2 PI - Load Change

image_folder = strcat(root_image_folder, '/Theorem_2_PI/Load_Change');
[~,~]=mkdir(image_folder);
image_folder = strcat(image_folder, '/');

theo2pi_config = default_config;
theo2pi_config.opt_theorem = 2;
theo2pi_config.opt_equilibrium_controller = 1;
theo2pi_config.image_folder = image_folder;
theo2pi_config.disturbance_Ro_enable = 1;
theo2pi_config.simulation_duration = 0.7;

% Buck
new_sim = theo2pi_config;
new_sim.circuit = circuit_buck;
new_sim.test_voltages = new_sim.circuit.single_voltage;
config_simulations{end+1} = new_sim;

% Boost
new_sim = theo2pi_config;
new_sim.circuit = circuit_boost;
new_sim.test_voltages = new_sim.circuit.single_voltage;
config_simulations{end+1} = new_sim;

% Buck-Boost
new_sim = theo2pi_config;
new_sim.circuit = circuit_buck_boost;
new_sim.test_voltages = new_sim.circuit.single_voltage;
config_simulations{end+1} = new_sim;


%% Run PWM - Input Change

image_folder = strcat(root_image_folder, '/PWM/Input_Change');
[~,~]=mkdir(image_folder);
image_folder = strcat(image_folder, '/');

pwm_config = default_config;
pwm_config.opt_pwm = 1;
pwm_config.image_folder = image_folder;
pwm_config.disturbance_Vin_enable = 1;
pwm_config.simulation_duration = 0.6;

% Buck
new_sim = pwm_config;
new_sim.circuit = circuit_buck;
new_sim.test_voltages = new_sim.circuit.single_voltage;
config_simulations{end+1} = new_sim;

% Boost
new_sim = pwm_config;
new_sim.circuit = circuit_boost;
new_sim.test_voltages = new_sim.circuit.single_voltage;
config_simulations{end+1} = new_sim;

% Buck-Boost
new_sim = pwm_config;
new_sim.circuit = circuit_buck_boost;
new_sim.test_voltages = new_sim.circuit.single_voltage;
config_simulations{end+1} = new_sim;

%% Run Theorem 1 - Input Change

image_folder = strcat(root_image_folder, '/Theorem_1/Input_Change');
[~,~]=mkdir(image_folder);
image_folder = strcat(image_folder, '/');

theo1_config = default_config;
theo1_config.opt_theorem = 1;
theo1_config.image_folder = image_folder;
theo1_config.disturbance_Vin_enable = 1;
theo1_config.disturbance_Vin_time = 0.25;
theo1_config.simulation_duration = 0.6;

% Buck
new_sim = theo1_config;
new_sim.circuit = circuit_buck;
new_sim.test_voltages = new_sim.circuit.single_voltage;
config_simulations{end+1} = new_sim;

% Boost
new_sim = theo1_config;
new_sim.circuit = circuit_boost;
new_sim.test_voltages = new_sim.circuit.single_voltage;
config_simulations{end+1} = new_sim;

% Buck-Boost
new_sim = theo1_config;
new_sim.circuit = circuit_buck_boost;
new_sim.test_voltages = new_sim.circuit.single_voltage;
config_simulations{end+1} = new_sim;


%% Run Theorem 2 - Input Change

image_folder = strcat(root_image_folder, '/Theorem_2/Input_Change');
[~,~]=mkdir(image_folder);
image_folder = strcat(image_folder, '/');

theo2_config = default_config;
theo2_config.opt_theorem = 2;
theo2_config.image_folder = image_folder;
theo2_config.disturbance_Vin_enable = 1;

% Buck
new_sim = theo2_config;
new_sim.circuit = circuit_buck;
new_sim.test_voltages = new_sim.circuit.single_voltage;
new_sim.disturbance_Vin_time = 0.2;
new_sim.simulation_duration = 0.45;
config_simulations{end+1} = new_sim;

% Boost
new_sim = theo2_config;
new_sim.circuit = circuit_boost;
new_sim.test_voltages = new_sim.circuit.single_voltage;
new_sim.disturbance_Vin_time = 1.8;
new_sim.simulation_duration = 2;
config_simulations{end+1} = new_sim;

% Buck-Boost
new_sim = theo2_config;
new_sim.circuit = circuit_buck_boost;
new_sim.test_voltages = new_sim.circuit.single_voltage;
new_sim.disturbance_Vin_time = 1.8;
new_sim.simulation_duration = 2;
config_simulations{end+1} = new_sim;


%% Run Theorem 2 PI - Input Change

image_folder = strcat(root_image_folder, '/Theorem_2_PI/Input_Change');
[~,~]=mkdir(image_folder);
image_folder = strcat(image_folder, '/');

theo2pi_config = default_config;
theo2pi_config.opt_theorem = 2;
theo2pi_config.opt_equilibrium_controller = 1;
theo2pi_config.image_folder = image_folder;
theo2pi_config.disturbance_Vin_enable = 1;
theo2pi_config.simulation_duration = 0.7;

% Buck
new_sim = theo2pi_config;
new_sim.circuit = circuit_buck;
new_sim.test_voltages = new_sim.circuit.single_voltage;
config_simulations{end+1} = new_sim;

% Boost
new_sim = theo2pi_config;
new_sim.circuit = circuit_boost;
new_sim.test_voltages = new_sim.circuit.single_voltage;
config_simulations{end+1} = new_sim;

% Buck-Boost
new_sim = theo2pi_config;
new_sim.circuit = circuit_buck_boost;
new_sim.test_voltages = new_sim.circuit.single_voltage;
config_simulations{end+1} = new_sim;

%% RUN EVERYTHING !!!!!!!


% Number of simulations to run
sim_N = length(config_simulations);

main_bar = waitbar(0, 'Prepare the computer!!', 'name', 'Running');

try
    for sim_id=1:sim_N

        sim_config = config_simulations{sim_id};

        opt_model = sim_config.opt_model;
        opt_discrete = sim_config.opt_discrete;
        opt_theorem = sim_config.opt_theorem;
        opt_pwm = sim_config.opt_pwm;
        opt_equilibrium_controller = sim_config.opt_equilibrium_controller;
        opt_update_equilibrium = sim_config.opt_update_equilibrium;
        opt_partial_information = sim_config.opt_partial_information;
        opt_constant_reference = sim_config.opt_constant_reference;
        disturbance_Vin_enable = sim_config.disturbance_Vin_enable;
        disturbance_Ro_enable = sim_config.disturbance_Ro_enable;
        disturbance_Vin_time = sim_config.disturbance_Vin_time;
        disturbance_Ro_time = sim_config.disturbance_Ro_time;
        Tsw = sim_config.switching_period;
        circuit = sim_config.circuit;
        test_voltages = sim_config.test_voltages;
        simulation_duration = sim_config.simulation_duration;
        image_folder = sim_config.image_folder;


        % Update wait bar
        waitbar((sim_id)/sim_N, main_bar, sprintf('Case %i of %d',sim_id, sim_N));

        run simulate_converter

        close all;

    end
catch exception
    close(main_bar);
    rethrow(exception);
end
close(main_bar);


