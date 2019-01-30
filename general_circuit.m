%% Initial Setup
clear; clc; close all;

addpath(genpath('functions'))
addpath(genpath('models'))

image_folder = 'images';

[~,~]=mkdir(image_folder);
image_folder = strcat(image_folder, '\');

plot_compression_rate = 1e3;

%% Desired Theorem to us

% Theorems defineds
%   1 - Fixed Equilibrium
%   2 - Valiable Equilibrium
opt_theorem = 2;

%% Desired DC-DC converter to use

% Options can be found in the system directory:
%   buck
%   boost
%   buck_boost
circuit = buck_boost;


R  = 2; % [Ohm]
Ro = 50; % [Ohm]
Co = 470e-6; % [F]
L  = 500e-6; % [H]

% Get the converter Space State representation and store it in the sys
% struct. The struct contains the fields
%   A - Cell Array with each subsystem A matrix - sys.A{i}
%   B - Cell Array with each subsystem B matrix - sys.B{i}
%   C - Cell Array with each subsystem C matrix - sys.C{i}
%   D - Cell Array with each subsystem D matrix - sys.D{i}
%   Q - Cell Array with each subsystem Q matrix - sys.Q{i}
%   Q - Cell Array with each subsystem Q matrix - sys.Q{i}
%   N - Number of avaliable subsystems - sys.N
sys = circuit.get_sys(R, Ro, Co, L);

%% System specifications

% Input voltage
Vs = 100; % [V]

% System starts discharged
x0 = [0; 0];

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
        % test_voltages = 100:10:200;
        % Buck
        % test_voltages = 10:10:90;
        % Buck-Boost
        test_voltages = 10:10:130;
        
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

    
    waitbar((Ns-i)/Ns, bar, sprintf('Simulation %i of %d',Ns-i+1, Ns));
    
    sim('general_system.slx', [0 0.1]);
    
    % Store only samples of the data, this will be made in order to save
    % memory use
    sim_out(i).x = compress_data(logsout.get('x').Values, plot_compression_rate);
end
close(bar);

plot_voltage_current(sim_out, circuit.name, image_folder);

plot_voltage_time(sim_out, circuit.name, image_folder);
