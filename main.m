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


plot_compression_rate = 1e3;


Simulink.fileGenControl('set', 'CacheFolder', cache_folder);

%% System Specifications

run system_specifications

%% Simulation Parametersr

% Model type to simulate
% Types of simulations:
%   1 - General Space State System Models
%   2 - Specific Circuit Model
opt_model = 2;

% Desired Theorem to use
% Theorems defines
%   1 - Fixed Equilibrium
%   2 - Valiable Equilibrium
opt_theorem = 2;

% Use PWM Controled mode or default switched control
% Options
%   0 - Use default control system
%   1 - Use pwm control system
opt_pwm = 0;


% Update the equilibrium point from the system
% Options
%   0 - Use given equilibrium
%   1 - Update equilibrium based on given reference voltage
opt_update_equilibrium = 1;

disturbance_Vin_enable = 0;
disturbance_Ro_enable = 0;

opt_constant_reference = 1;

% Desired DC-DC converter to use
% Options can be found in the system directory:
%   buck
%   boost
%   buck_boost
%   buck_boost_non_inverting
circuit = buck_boost_non_inverting(R, Ro, Co, L);


test_voltages = circuit.test_voltages;

test_voltages = [300];

simulation_duration = 1;


%% Prepare Data

pwm_pid_kp = circuit.pwm_pid_kp;
pwm_pid_ki = circuit.pwm_pid_ki;
pwm_pid_kd = circuit.pwm_pid_kd;

reference_pid_kp = circuit.reference_pid_kp;
reference_pid_ki = circuit.reference_pid_ki;

[reference_ve_limit_lower, reference_ve_limit_upper] = get_reference_ve_limits(circuit, Vs);

run load_circuit_sys

run circuit_disturbance

%% Lambdas to simulate

% lambdas = generate_lambda_voltage(sys, test_voltages);

step = 0.25;
sequence = 0:step:1;
test_lambdas = zeros(15, 3);
index = 1;
for lambda1=sequence
    for lambda2=0:step:(1-lambda1)
        lambda3 = 1 - lambda2 - lambda1;
        test_lambdas(index,1) = lambda1;
        test_lambdas(index,2) = lambda2;
        test_lambdas(index,3) = lambda3;
        index = index+1;
    end
end
lambdas = test_lambdas;

lambdas = [0.1 0.2 0.7];


% lambdas = generate_lambda_voltage(sys, test_voltages);

%% Simulate Converter 

% Get model to simulate
switch(opt_model)
    case 1
        model = 'general_system.slx';
    case 2
        model = circuit.simulink;
end

% Number of simulations to run
Ns = size(lambdas, 1);

bar = waitbar(0, 'Preparing simulation', 'name', 'Simulating');

load_system(model);

for i=Ns:-1:1
    
    % Convert the truct used to represent the space state to double, so it
    % can be used in the simulink
    [A, B, C, D, Q] = gss2double(sys);
    
    % Calculate the P matrix, as the equilibrium point. Calculation will be
    % in accordance with the chosen control theorem
    switch(opt_theorem)
        case 1
            [P, xe] = calc_sys_theorem_1(sys, lambdas(i,:));
        case 2
            [P, xe] = calc_sys_theorem_2(sys, lambdas(i,:));
    end
    
    xe = [0 40]';
    
    % Creates a bus, wich will be used in the simulink to simplify the
    % model
    SystemDataBus = create_bus_SystemDataBus(A, B, P, Q, sys.N);

    % Update wait bar
    waitbar((Ns-i)/Ns, bar, sprintf('Simulation %i of %d',Ns-i+1, Ns));
    
    % Run simulation
    sim(model, simulation_duration);
    
    % Store only samples of the data, this will be made in order to save
    % memory use
    sim_out(i).Iout = compress_data(logsout.get('IL').Values, plot_compression_rate);
    sim_out(i).Vout = compress_data(logsout.get('Vout').Values, plot_compression_rate);
    sim_out(i).xe = compress_data(logsout.get('xe').Values, plot_compression_rate);
end
close(bar);

%% Analysis

plot_voltage_current(sim_out, circuit.name, image_folder);

plot_voltage_time(sim_out, circuit.name, image_folder);

plot_current_time(sim_out, circuit.name, image_folder);

if disturbance_Ro_enable || disturbance_Vin_enable
    plot_disturbance_voltage_time(sim_out, disturbance_Ro_time, circuit.name, image_folder);
end