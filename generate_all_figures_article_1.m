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
addpath(genpath('data'))


plot_compression_rate = 1e0;


Simulink.fileGenControl('set', 'CacheFolder', cache_folder);
Simulink.fileGenControl('set', 'CodeGenFolder', cache_folder);

warning('off','Simulink:Engine:UINotUpdatedDuringRapidAccelSim')


image_folder = strcat(root_image_folder, '/Article');
[~,~]=mkdir(image_folder);
image_folder = strcat(image_folder, '/');

root_data_mat_folder = 'data_mat';

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
saveVars = {};

default_config = simulation_config();
default_config.opt_model = 2;
default_config.opt_discrete = false;
default_config.opt_theorem = 2;
default_config.opt_range_design = true;
default_config.opt_pwm = false;
default_config.opt_pwm_current_controller = false;
default_config.opt_update_equilibrium = true;
default_config.opt_equilibrium_controller = false;
default_config.opt_current_correction = false;
default_config.opt_partial_information = false;
default_config.opt_constant_reference = true;
default_config.variable_reference_step_period = 0.5;
default_config.disturbance_Vin_enable = false;
default_config.disturbance_Ro_enable = false;
default_config.disturbance_Vin_time = disturbance_Vin_time;
default_config.disturbance_Ro_time = disturbance_Ro_time;
default_config.opt_variable_load = false;
default_config.opt_dead_time = true;
default_config.opt_mode_hopping = false;
default_config.opt_sensor_noises = false;
default_config.switching_period = 1/40e3;
default_config.simulation_sample = 1e-7;
default_config.circuit = circuit_buck_boost;
default_config.test_voltages = circuit_buck_boost.single_voltage;
default_config.simulation_duration = 0.15;
default_config.opt_measurement_frequency = false;
default_config.opt_measurement_efficiency = false;
default_config.opt_measurement_clock = false;
default_config.image_folder = image_folder;

opt_theorem_discrete = 3;

%% Run - Frequency Test - Ideal 1MHz - Continuous Controller 1

new_sim = copy(default_config);
new_sim.opt_theorem = 1;
new_sim.opt_measurement_frequency = true;
new_sim.opt_measurement_ripple = true;
new_sim.opt_measurement_error = true;
new_sim.opt_constant_reference = false;
new_sim.opt_dead_time = false;
new_sim.switching_period = 1/1000e3;
new_sim.simulation_sample = 1e-7;
new_sim.simulation_duration = 12.5;

config_simulations{end+1} = new_sim;
continuous_1_frequency_ideal_1M = new_sim;
saveVars.continuous_1_frequency_ideal_1M = new_sim;


%% Run - Frequency Test - Ideal 1MHz - Continuous Controller 2

new_sim = copy(default_config);
new_sim.opt_measurement_frequency = true;
new_sim.opt_measurement_ripple = true;
new_sim.opt_measurement_error = true;
new_sim.opt_constant_reference = false;
new_sim.variable_reference_step_period = 0.5;
new_sim.opt_dead_time = false;
new_sim.switching_period = 1/1000e3;
new_sim.simulation_sample = 1e-7;
new_sim.simulation_duration = 12.5;

config_simulations{end+1} = new_sim;
continuous_2_frequency_ideal_1M = new_sim;
saveVars.continuous_2_frequency_ideal_1M = new_sim;


%% Run - Frequency Test - Ideal 1MHz - Discrete Controller

new_sim = copy(default_config);
new_sim.opt_discrete = true;
new_sim.opt_theorem = opt_theorem_discrete;
new_sim.opt_measurement_frequency = true;
new_sim.opt_measurement_ripple = true;
new_sim.opt_measurement_error = true;
new_sim.opt_constant_reference = false;
new_sim.opt_dead_time = false;
new_sim.switching_period = 1/1000e3;
new_sim.simulation_sample = 1e-7;
new_sim.simulation_duration = 12.5;

config_simulations{end+1} = new_sim;
discrete_frequency_ideal_1M = new_sim;
saveVars.discrete_frequency_ideal_1M = new_sim;


%% Run - Frequency Test - Ideal 200kHz - Continuous Controller 1

new_sim = copy(default_config);
new_sim.opt_theorem = 1;
new_sim.opt_measurement_frequency = true;
new_sim.opt_measurement_ripple = true;
new_sim.opt_measurement_error = true;
new_sim.opt_constant_reference = false;
new_sim.opt_dead_time = false;
new_sim.switching_period = 1/200e3;
new_sim.simulation_sample = 1e-7;
new_sim.simulation_duration = 12.5;

config_simulations{end+1} = new_sim;
continuous_1_frequency_ideal_200 = new_sim;
saveVars.continuous_1_frequency_ideal_200 = new_sim;


%% Run - Frequency Test - Ideal 200kHz - Continuous Controller 2

new_sim = copy(default_config);
new_sim.opt_measurement_frequency = true;
new_sim.opt_measurement_ripple = true;
new_sim.opt_measurement_error = true;
new_sim.opt_constant_reference = false;
new_sim.variable_reference_step_period = 0.5;
new_sim.opt_dead_time = false;
new_sim.switching_period = 1/200e3;
new_sim.simulation_sample = 1e-7;
new_sim.simulation_duration = 12.5;

config_simulations{end+1} = new_sim;
continuous_2_frequency_ideal_200 = new_sim;
saveVars.continuous_2_frequency_ideal_200 = new_sim;


%% Run - Frequency Test - Ideal 200kHz - Discrete Controller

new_sim = copy(default_config);
new_sim.opt_discrete = true;
new_sim.opt_theorem = opt_theorem_discrete;
new_sim.opt_measurement_frequency = true;
new_sim.opt_measurement_ripple = true;
new_sim.opt_measurement_error = true;
new_sim.opt_constant_reference = false;
new_sim.opt_dead_time = false;
new_sim.switching_period = 1/200e3;
new_sim.simulation_sample = 1e-7;
new_sim.simulation_duration = 12.5;

config_simulations{end+1} = new_sim;
discrete_frequency_ideal_200 = new_sim;
saveVars.discrete_frequency_ideal_200 = new_sim;


%% Run - Frequency Test - Ideal 40kHz - Continuous Controller 1

new_sim = copy(default_config);
new_sim.opt_theorem = 1;
new_sim.opt_measurement_frequency = true;
new_sim.opt_measurement_ripple = true;
new_sim.opt_measurement_error = true;
new_sim.opt_constant_reference = false;
new_sim.opt_dead_time = false;
new_sim.switching_period = 1/40e3;
new_sim.simulation_duration = 12.5;

config_simulations{end+1} = new_sim;
continuous_1_frequency_ideal_40 = new_sim;
saveVars.continuous_1_frequency_ideal_40 = new_sim;


%% Run - Frequency Test - Ideal 40kHz - Continuous Controller 2

new_sim = copy(default_config);
new_sim.opt_measurement_frequency = true;
new_sim.opt_measurement_ripple = true;
new_sim.opt_measurement_error = true;
new_sim.opt_constant_reference = false;
new_sim.variable_reference_step_period = 0.5;
new_sim.opt_dead_time = false;
new_sim.switching_period = 1/40e3;
new_sim.simulation_duration = 12.5;

config_simulations{end+1} = new_sim;
continuous_2_frequency_ideal_40 = new_sim;
saveVars.continuous_2_frequency_ideal_40 = new_sim;


%% Run - Frequency Test - Ideal 40kHz - Discrete Controller

new_sim = copy(default_config);
new_sim.opt_discrete = true;
new_sim.opt_theorem = opt_theorem_discrete;
new_sim.opt_measurement_frequency = true;
new_sim.opt_measurement_ripple = true;
new_sim.opt_measurement_error = true;
new_sim.opt_constant_reference = false;
new_sim.opt_dead_time = false;
new_sim.switching_period = 1/40e3;
new_sim.simulation_duration = 12.5;

config_simulations{end+1} = new_sim;
discrete_frequency_ideal_40 = new_sim;
saveVars.discrete_frequency_ideal_40 = new_sim;


%% Run - Frequency Test - Ideal 40kHz - PI - Continuous Controller 1

new_sim = copy(default_config);
new_sim.opt_theorem = 1;
new_sim.opt_current_correction = true;
new_sim.opt_measurement_frequency = true;
new_sim.opt_measurement_ripple = true;
new_sim.opt_measurement_error = true;
new_sim.opt_constant_reference = false;
new_sim.opt_dead_time = false;
new_sim.switching_period = 1/40e3;
new_sim.simulation_duration = 12.5;

config_simulations{end+1} = new_sim;
continuous_1_frequency_ideal_PI_40 = new_sim;
saveVars.continuous_1_frequency_ideal_PI_40 = new_sim;


%% Run - Frequency Test - Ideal 40kHz - PI - Continuous Controller 2

new_sim = copy(default_config);
new_sim.opt_current_correction = true;
new_sim.opt_measurement_frequency = true;
new_sim.opt_measurement_ripple = true;
new_sim.opt_measurement_error = true;
new_sim.opt_constant_reference = false;
new_sim.opt_dead_time = false;
new_sim.switching_period = 1/40e3;
new_sim.simulation_duration = 12.5;

config_simulations{end+1} = new_sim;
continuous_2_frequency_ideal_PI_40 = new_sim;
saveVars.continuous_2_frequency_ideal_PI_40 = new_sim;


%% Run - Frequency Test - Ideal 40kHz - PI - Discrete Controller

new_sim = copy(default_config);
new_sim.opt_discrete = true;
new_sim.opt_theorem = opt_theorem_discrete;
new_sim.opt_current_correction = true;
new_sim.opt_measurement_frequency = true;
new_sim.opt_measurement_ripple = true;
new_sim.opt_measurement_error = true;
new_sim.opt_constant_reference = false;
new_sim.opt_dead_time = false;
new_sim.switching_period = 1/40e3;
new_sim.simulation_duration = 12.5;

config_simulations{end+1} = new_sim;
discrete_frequency_ideal_PI_40 = new_sim;
saveVars.discrete_frequency_ideal_PI_40 = new_sim;


%% Run - Frequency Test - Real 40kHz - PI - Continuous Controller 1

new_sim = copy(default_config);
new_sim.opt_theorem = 1;
new_sim.opt_current_correction = true;
new_sim.opt_measurement_frequency = true;
new_sim.opt_measurement_ripple = true;
new_sim.opt_measurement_error = true;
new_sim.opt_constant_reference = false;
new_sim.simulation_duration = 12.5;

config_simulations{end+1} = new_sim;
continuous_1_frequency_real_PI = new_sim;
saveVars.continuous_1_frequency_real_PI = new_sim;


%% Run - Frequency Test - Real 40kHz - PI - Continuous Controller 2

new_sim = copy(default_config);
new_sim.opt_current_correction = true;
new_sim.opt_measurement_frequency = true;
new_sim.opt_measurement_ripple = true;
new_sim.opt_measurement_error = true;
new_sim.opt_constant_reference = false;
new_sim.simulation_duration = 12.5;

config_simulations{end+1} = new_sim;
continuous_2_frequency_real_PI = new_sim;
saveVars.continuous_2_frequency_real_PI = new_sim;


%% Run 5 - Frequency Test - Real 40kHz - PI - Discrete Controller

new_sim = copy(default_config);
new_sim.opt_discrete = true;
new_sim.opt_theorem = opt_theorem_discrete;
new_sim.opt_current_correction = true;
new_sim.opt_measurement_frequency = true;
new_sim.opt_measurement_ripple = true;
new_sim.opt_measurement_error = true;
new_sim.opt_constant_reference = false;
new_sim.simulation_duration = 12.5;

config_simulations{end+1} = new_sim;
discrete_frequency_real_PI = new_sim;
saveVars.discrete_frequency_real_PI = new_sim;


%% Run - Load Step - Classic Controller

new_sim = copy(default_config);
new_sim.opt_pwm = true;
new_sim.disturbance_Ro_enable = true;
new_sim.disturbance_Ro_time = 2;
new_sim.simulation_duration = 3.5;

config_simulations{end+1} = new_sim;
classic_load_step = new_sim;
saveVars.classic_load_step = new_sim;


%% Run - Load Step - Continuous Controller 1

new_sim = copy(default_config);
new_sim.opt_theorem = 1;
new_sim.opt_current_correction = true;
new_sim.disturbance_Ro_enable = true;
new_sim.disturbance_Ro_time = 2;
new_sim.test_voltages = new_sim.circuit.single_voltage;
new_sim.simulation_duration = 3.5;

config_simulations{end+1} = new_sim;
continuous_1_load_step = new_sim;
saveVars.continuous_1_load_step = new_sim;


%% Run - Load Step - Continuous Controller 2

new_sim = copy(default_config);
new_sim.opt_current_correction = true;
new_sim.disturbance_Ro_enable = true;
new_sim.disturbance_Ro_time = 2;
new_sim.test_voltages = new_sim.circuit.single_voltage;
new_sim.simulation_duration = 3.5;

config_simulations{end+1} = new_sim;
continuous_2_load_step = new_sim;
saveVars.continuous_2_load_step = new_sim;

%% Run - Load Step - Discrete Controller

new_sim = copy(default_config);
new_sim.opt_discrete = true;
new_sim.opt_theorem = opt_theorem_discrete;
new_sim.opt_current_correction = true;
new_sim.disturbance_Ro_enable = true;
new_sim.disturbance_Ro_time = 2;
new_sim.test_voltages = new_sim.circuit.single_voltage;
new_sim.simulation_duration = 3.5;

config_simulations{end+1} = new_sim;
discrete_load_step = new_sim;
saveVars.discrete_load_step = new_sim;

%% Run - Reference Step - Continuous Controller 1

new_sim = copy(default_config);
new_sim.opt_theorem = 1;
new_sim.test_voltages = new_sim.circuit.single_voltage;
new_sim.simulation_duration = 1;

config_simulations{end+1} = new_sim;
continuous_1_reference_step = new_sim;
saveVars.continuous_1_reference_step = new_sim;

%% Run - Reference Step - Continuous Controller 2

new_sim = copy(default_config);
new_sim.test_voltages = new_sim.circuit.single_voltage;
new_sim.simulation_duration = 1;

config_simulations{end+1} = new_sim;
continuous_2_reference_step = new_sim;
saveVars.continuous_2_reference_step = new_sim;


%% Run - Reference Step - Discrete Controller

new_sim = copy(default_config);
new_sim.opt_discrete = true;
new_sim.opt_theorem = opt_theorem_discrete;
new_sim.test_voltages = new_sim.circuit.single_voltage;
new_sim.simulation_duration = 1;

config_simulations{end+1} = new_sim;
discrete_reference_step = new_sim;
saveVars.discrete_reference_step = new_sim;


%% Run - Reference Step - Classic Controller

new_sim = copy(default_config);
new_sim.opt_pwm = true;
new_sim.test_voltages = new_sim.circuit.single_voltage;
new_sim.simulation_duration = 1;

config_simulations{end+1} = new_sim;
classic_reference_step = new_sim;
saveVars.classic_reference_step = new_sim;


%% Run 11 - Ripple Test - Classic Controller

new_sim = copy(default_config);
new_sim.opt_pwm = true;
new_sim.opt_measurement_ripple = true;
new_sim.opt_measurement_error = true;
new_sim.opt_constant_reference = false;
new_sim.test_voltages = new_sim.circuit.single_voltage;
new_sim.simulation_duration = 12.5;

config_simulations{end+1} = new_sim;
classic_ripple = new_sim;
saveVars.classic_ripple = new_sim;

%% Run - Step Analysis - Continuous Controller 1

new_sim = copy(default_config);
new_sim.opt_theorem = 1;
new_sim.test_voltages = new_sim.circuit.operation_range_voltage_min:5:new_sim.circuit.operation_range_voltage_max;
new_sim.simulation_duration = 0.5;

config_simulations{end+1} = new_sim;
continuous_1_step_analysis = new_sim;
saveVars.continuous_1_step_analysis = new_sim;

%% Run - Step Analysis - Continuous Controller 2

new_sim = copy(default_config);
new_sim.opt_theorem = 2;
new_sim.test_voltages = new_sim.circuit.operation_range_voltage_min:5:new_sim.circuit.operation_range_voltage_max;
new_sim.simulation_duration = 1;

config_simulations{end+1} = new_sim;
continuous_2_step_analysis = new_sim;
saveVars.continuous_2_step_analysis = new_sim;


%% Run - Step Analysis - Discrete Controller

new_sim = copy(default_config);
new_sim.opt_discrete = true;
new_sim.opt_theorem = opt_theorem_discrete;
new_sim.test_voltages = new_sim.circuit.operation_range_voltage_min:5:new_sim.circuit.operation_range_voltage_max;
new_sim.simulation_duration = 0.5;

config_simulations{end+1} = new_sim;
discrete_step_analysis = new_sim;
saveVars.discrete_step_analysis = new_sim;


%% Run - Step Analysis - Classic Controller

new_sim = copy(default_config);
new_sim.opt_pwm = true;
new_sim.test_voltages = new_sim.circuit.operation_range_voltage_min:5:new_sim.circuit.operation_range_voltage_max;
new_sim.simulation_duration = 0.5;

config_simulations{end+1} = new_sim;
classic_step_analysis = new_sim;
saveVars.classic_step_analysis = new_sim;


%% Run - Step Analysis - Continuous Controller 1 - Single Design

new_sim = copy(default_config);
new_sim.opt_theorem = 1;
new_sim.opt_range_design = false;
new_sim.test_voltages = new_sim.circuit.operation_range_voltage_min:5:new_sim.circuit.operation_range_voltage_max;
new_sim.simulation_duration = 0.5;

config_simulations{end+1} = new_sim;
continuous_1_step_analysis_single = new_sim;
saveVars.continuous_1_step_analysis_single = new_sim;

%% Run - Step Analysis - Continuous Controller 2 - Single Design

new_sim = copy(default_config);
new_sim.opt_theorem = 2;
new_sim.opt_range_design = false;
new_sim.test_voltages = new_sim.circuit.operation_range_voltage_min:5:new_sim.circuit.operation_range_voltage_max;
new_sim.simulation_duration = 1;

config_simulations{end+1} = new_sim;
continuous_2_step_analysis_single = new_sim;
saveVars.continuous_2_step_analysis_single = new_sim;


%% Run - Step Analysis - Discrete Controller - Single Design

new_sim = copy(default_config);
new_sim.opt_discrete = true;
new_sim.opt_theorem = opt_theorem_discrete;
new_sim.opt_range_design = false;
new_sim.test_voltages = new_sim.circuit.operation_range_voltage_min:5:new_sim.circuit.operation_range_voltage_max;
new_sim.simulation_duration = 0.5;

config_simulations{end+1} = new_sim;
discrete_step_analysis_single = new_sim;
saveVars.discrete_step_analysis_single = new_sim;


%% RUN EVERYTHING !!!!!!!


% Number of simulations to run
sim_N = length(config_simulations);

main_bar = waitbar(0, 'Prepare the computer!!', 'name', 'Running');

total_time = 0;
sim_time = 0;
exec_time = 0;
for i=1:sim_N
    sim_config = config_simulations{i};
    Ns = length(sim_config.test_voltages);
    total_time = total_time + Ns*sim_config.simulation_duration;
end

try
    for sim_id=1:sim_N
        tic
        sim_config = config_simulations{sim_id};

        opt_model = sim_config.opt_model;
        opt_discrete = sim_config.opt_discrete;
        opt_theorem = sim_config.opt_theorem;
        opt_range_design = sim_config.opt_range_design;
        opt_pwm = sim_config.opt_pwm;
        opt_pwm_current_controller = sim_config.opt_pwm_current_controller;
        opt_update_equilibrium = sim_config.opt_update_equilibrium;
        opt_equilibrium_controller = sim_config.opt_equilibrium_controller;
        opt_current_correction = sim_config.opt_current_correction;
        opt_partial_information = sim_config.opt_partial_information;
        opt_constant_reference = sim_config.opt_constant_reference;
        variable_reference_step_period = sim_config.variable_reference_step_period;
        disturbance_Vin_enable = sim_config.disturbance_Vin_enable;
        disturbance_Ro_enable = sim_config.disturbance_Ro_enable;
        disturbance_Vin_time = sim_config.disturbance_Vin_time;
        disturbance_Ro_time = sim_config.disturbance_Ro_time;
        opt_variable_load = sim_config.opt_variable_load;
        opt_dead_time = sim_config.opt_dead_time;
        opt_mode_hopping = sim_config.opt_mode_hopping;
        opt_sensor_noises = sim_config.opt_sensor_noises;
        Ts = sim_config.switching_period;
        Ti = sim_config.simulation_sample;
        mu = -1;
        circuit = sim_config.circuit;
        test_voltages = sim_config.test_voltages;
        simulation_duration = sim_config.simulation_duration;
        opt_measurement_frequency = sim_config.opt_measurement_frequency;
        opt_measurement_efficiency = sim_config.opt_measurement_efficiency;
        opt_measurement_ripple = sim_config.opt_measurement_ripple;
        opt_measurement_clock = sim_config.opt_measurement_clock;
        opt_measurement_error = sim_config.opt_measurement_error;
        image_folder = sim_config.image_folder;


        % Update wait bar
        if isvalid(main_bar)
            if sim_time>0
                val = (total_time - sim_time) * exec_time/sim_time;
%                 val = (sim_N - sim_id) * exec_time/(sim_id-1);
                if val < 60
                    remaining = sprintf('%.fs', val);
                else
                    remaining = sprintf('%.fm %.fs', floor(val/60), mod(val,60));
                end
            else
                remaining = 'unknown';
            end
            text = sprintf('Case %i of %d - Simulated %.fs of %.fs\nTime remaining: %s',sim_id, sim_N, sim_time, total_time, remaining);
            waitbar(sim_time/total_time, main_bar, text);
%             waitbar((sim_id-1)/sim_N, main_bar, text);
        end
        
        run simulate_converter
        close all;
        
        config_simulations{sim_id}.sim_out = sim_out;
        
        Ns = length(sim_config.test_voltages);
        sim_time = sim_time + Ns*sim_config.simulation_duration;
        exec_time = exec_time + toc;
    end
catch exception
    close(main_bar);
    rethrow(exception);
end
close(main_bar);

%% End of simulations

file = fullfile(root_data_mat_folder, 'all_simulations_article.mat');

if isfile(file)
    save(file, '-struct', 'saveVars')
else
    save(file, '-struct', 'saveVars')
end

%%

run generate_all_figures_article_2

load handel
sound(y,Fs)
