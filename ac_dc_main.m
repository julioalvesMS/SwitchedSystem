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

run ac_dc_system_specifications

%% Simulation Parametersr

% Model type to simulate
% Types of simulations:
%   1 - General Space State System Models
%   2 - Specific Circuit Model
opt_model = 2;

% Choose between continuous or discrete model
% Options:
%   0 - Continuous Controller
%   1 - Discrete Controller
opt_discrete = false;

% Desired Theorem to use
% Theorems defines
%   1 - Fixed Equilibrium
%   2 - Valiable Equilibrium
opt_theorem = 2;

% Use PWM Controled mode or default switched control
% Options
%   0 - Use default control system
%   1 - Use pwm control system
opt_pwm = false;

% Use Current Control or only Voltage Control for PWM
% Options
%   0 - Use only voltage control
%   1 - Use voltage and current control
opt_pwm_current_controller = false;

% Choose between a constant output voltage or one with a different profile
% Options
%   0 - Update reference according to the profile in the simulink
%   1 - Use constante reference
opt_constant_reference = true;

% Disturbances to be applied during simulations
% Options
%   disturbance_Vin_enable - Enable step disturbance in the input voltage
%   disturbance_Ro_enable - Enable step disturbance in the load resistance
disturbance_Vin_enable = false;
disturbance_Ro_enable = false;

% Simulate the dead time observed in real life switches
% Options
%   0 - Consider ideal switches
%   1 - Consider dead time
opt_dead_time = false;

% Desired AC-DC converter to use
% Options can be found in the system directory:
%   rectifier
circuit = rectifier(R, L, Co);

test_voltages = circuit.single_voltage;


simulation_duration = 0.4;


%% Prepare Data

pwm_pid_kp = circuit.pwm_pid_kp;
pwm_pid_ki = circuit.pwm_pid_ki;
        
pwm_pid_vc_vp = circuit.pwm_pid_vc_vp;
pwm_pid_vc_vi = circuit.pwm_pid_vc_vi;
pwm_pid_vc_cp = circuit.pwm_pid_vc_cp;
pwm_pid_vc_ci = circuit.pwm_pid_vc_ci;

run ac_dc_load_circuit_sys

%% Lambdas to simulate

lambdas = ones(1,sys.N)*(1/sys.N);

%% Simulate Converter 

% Get model to simulate
switch(opt_model)
    case 1
        if opt_discrete == 0
            model = 'general_system.slx';
        else
            model = 'discrete_general_system.slx';
        end
    case 2
        if opt_discrete == 0
            model = 'sim_rectifier.slx';
        else
            model = 'sim_discrete.slx';
        end
end

% Number of simulations to run
Ns = size(lambdas, 1);

bar = waitbar(0, 'Preparing simulation', 'name', 'Simulating');

try
    load_system(model);

    for i=Ns:-1:1
        
        Vref = test_voltages(i);

        % Calculate the P matrix, as the equilibrium point. Calculation will be
        % in accordance with the chosen control theorem
        if opt_discrete == 0
            switch(opt_theorem)
                case 1
                    [P, xe] = calc_sys_theorem_1(sys, lambdas(i,:));
                case 2
                    [P, xe] = calc_sys_theorem_2(sys, lambdas(i,:), [Vref 0 0]');
            end
        else
            [P, h, d, xe, dsys] = calc_sys_discrete_theorem_1(sys, dsys, lambdas(i,:));
        end

        % Convert the truct used to represent the space state to double, so it
        % can be used in the simulink
        [A, B, C, D, Q] = gss2double(sys);          % Continuos
        [Ad, Bd, Cd, Dd, Qd, Ld] = gss2double(dsys);    % Discrete


        % Creates a bus, wich will be used in the simulink to simplify the
        % model
        if opt_discrete == 0
            SystemDataBus = create_bus_SystemDataBus(A, B, P, Q, sys.N);
        else
            SystemDataBus = create_bus_DiscreteSystemDataBus(Ad, Bd, P, Qd, sys.N, h, Ld);
        end
            

        % Update wait bar
        waitbar((Ns-i)/Ns, bar, sprintf('Simulation %i of %d',Ns-i+1, Ns));

        % Run simulation
        sim(model, simulation_duration);

        % Store only samples of the data, this will be made in order to save
        % memory use
        sim_out(i).Vabc = downsample_timeseries(logsout.get('Vabc').Values, plot_compression_rate);
        sim_out(i).Vdc = downsample_timeseries(logsout.get('Vdc').Values, plot_compression_rate);
        sim_out(i).Vref = downsample_timeseries(logsout.get('Vref').Values, plot_compression_rate);
    end
    
    
catch exception
    close(bar);
    rethrow(exception);
end
close(bar);

%% Analysis
