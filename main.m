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


plot_compression_rate = 1e3;


Simulink.fileGenControl('set', 'CacheFolder', cache_folder);

%% Simulation Parametersr

% Model type to simulate
% Types of simulations:
%   1 - General Space State System Models
%   2 - Specific Circuit Model
opt_model = 1;

% Desired Theorem to use
% Theorems defines
%   1 - Fixed Equilibrium
%   2 - Valiable Equilibrium
opt_theorem = 2;

% Use PWM Controled mode or default switched control
% Options
%   0 - Use default control system
%   1 - Use pwm control system
opt_pwm = 1;

% Desired DC-DC converter to use
% Options can be found in the system directory:
%   buck
%   boost
%   buck_boost
circuit = boost;

%% System Specifications

% Simulation specifications
Ti=1e-6; % [s] - Simulation pace
Fsw = 10e3; % [H]

% Circuit specifications
R  = 2; % [Ohm] - Converter Resistance
Ro = 50; % [Ohm] - Load Resistance
Co = 470e-6; % [F] - Load Capacitance
L  = 500e-6; % [H]

% Input voltage
Vs = 100; % [V]

% System starts discharged
x0 = [0; 0];

%% Prepare Data

% Get the converter Space State representation and store it in the sys
% struct. The struct contains the fields
%   A - Cell Array with each subsystem A matrix - sys.A{i}
%   B - Cell Array with each subsystem B matrix - sys.B{i}
%   C - Cell Array with each subsystem C matrix - sys.C{i}
%   D - Cell Array with each subsystem D matrix - sys.D{i}
%   Q - Cell Array with each subsystem Q matrix - sys.Q{i}
%   N - Number of avaliable subsystems - sys.N
sys = circuit.get_sys(R, Ro, Co, L);

% Complete the sys struct with this data
sys.U = Vs;
sys.x0 = x0;

%% Lambdas to simulate

switch(opt_theorem)
    case 1
        % Generates a generic lambda matrix to be used
        lambdas = generate_lambda_2d(0.05);
    case 2
        % Boost
        test_voltages = 100:10:200;
        % Buck
        %test_voltages = 10:10:90;
        % Buck-Boost
        %test_voltages = 10:10:130;
        
        [sample_lambdas, equilibrium] = generate_sample_points(sys);
        
        % To be able to interpolate the data, we will use the only the
        % first lambda and the output voltage
        lamb1 = sample_lambdas(:,1);
        Vlamb = equilibrium(:,2);
        
        ns = length(test_voltages);
        lambdas = zeros(ns, 2);
        for i=1:ns
            % Only use the lambdas before the inflection point from the
            % voltage
            [~, inflection] = max(Vlamb);
            
            % interpolate the data to discover the needed lambda to reach
            % the desired voltage
            lamb = interp1(Vlamb(1:inflection), lamb1(1:inflection), test_voltages(i));
            lambdas(i,:) = [lamb, 1-lamb];
        end
end

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
    
    % Creates a bus, wich will be used in the simulink to simplify the
    % model
    SystemDataBus = create_bus_SystemDataBus(A, B, P, Q, xe, sys.N);

    % Update wait bar
    waitbar((Ns-i)/Ns, bar, sprintf('Simulation %i of %d',Ns-i+1, Ns));
    
    % Run simulation
    sim(model, [0 0.2]);
    
    % Store only samples of the data, this will be made in order to save
    % memory use
    sim_out(i).Iout = compress_data(logsout.get('Iout').Values, plot_compression_rate);
    sim_out(i).Vout = compress_data(logsout.get('Vout').Values, plot_compression_rate);
end
close(bar);

%% Analysis

plot_voltage_current(sim_out, circuit.name, image_folder);

plot_voltage_time(sim_out, circuit.name, image_folder);

plot_current_time(sim_out, circuit.name, image_folder);