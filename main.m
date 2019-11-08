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


plot_compression_rate = 1e1;


Simulink.fileGenControl('set', 'CacheFolder', cache_folder);


s = tf('s')
%% System Specifications

run system_specifications

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
opt_discrete = 0;


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


% Update the equilibrium point from the system
% Options
%   0 - Use given equilibrium
%   1 - Update equilibrium based on given reference voltage
opt_update_equilibrium = 1;


% Use a PI to determine and update the equilibrium point
% Options
%   0 - Don't use the PI
%   1 - Update the voltage from the equilibrium point using a PI controller
opt_equilibrium_controller = 1;


opt_partial_information = 0;

% Choose between a constant output voltage or one with a different profile
% Options
%   0 - Update reference according to the profile in the simulink
%   1 - Use constante reference
opt_constant_reference = 1;


% Disturbances to be applied during simulations
% Options
%   disturbance_Vin_enable - Enable step disturbance in the input voltage
%   disturbance_Ro_enable - Enable step disturbance in the load resistance
disturbance_Vin_enable = 0;
disturbance_Ro_enable = 0;

% Desired DC-DC converter to use
% Options can be found in the system directory:
%   buck
%   boost
%   buck_boost
%   buck_boost_non_inverting
circuit = buck_boost(R, Ro, Co, L);


test_voltages = circuit.test_voltages;
test_voltages = circuit.single_voltage;

% test_voltages = [190];

simulation_duration = 0.5;


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

if length(sys.A) == 2
    lambdas = generate_lambda_voltage(sys, test_voltages);
else

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

    lambdas = [0.4 0.2 0.4];
end

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
            model = circuit.simulink;
        else
            model = circuit.discrete_simulink;
        end
end

% Number of simulations to run
Ns = size(lambdas, 1);

bar = waitbar(0, 'Preparing simulation', 'name', 'Simulating');

try
    load_system(model);

    for i=Ns:-1:1

        % Calculate the P matrix, as the equilibrium point. Calculation will be
        % in accordance with the chosen control theorem
        if opt_discrete == 0
            switch(opt_theorem)
                case 1
                    [P, xe] = calc_sys_theorem_1(sys, lambdas(i,:));
                case 2
                    [P, xe] = calc_sys_theorem_2(sys, lambdas(i,:));
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
        sim_out(i).IL = downsample(logsout.get('IL').Values, plot_compression_rate);
        sim_out(i).Vout = downsample(logsout.get('Vout').Values, plot_compression_rate);
        sim_out(i).xe = downsample(logsout.get('xe').Values, plot_compression_rate);
        sim_out(i).Vref = downsample(logsout.get('Vref').Values, plot_compression_rate);
    end
catch exception
    close(bar);
    rethrow(exception);
end
close(bar);

%% Analysis

plot_voltage_current(sim_out, circuit.name, image_folder);

plot_voltage_time(sim_out, circuit.name, image_folder);

plot_current_time(sim_out, circuit.name, image_folder);

if disturbance_Vin_enable == 1
    plot_disturbance_voltage_time(sim_out, disturbance_Ro_time, circuit.name, image_folder);
end

if disturbance_Ro_enable == 1
    plot_disturbance_voltage_time(sim_out, disturbance_Vin_time, circuit.name, image_folder);
end

% figure;
% switches = logsout.get('State').Values.Data;
% Fa = abs(fft(switches));
% F = Fa(1:round(end/2));
% f = 0:1/simulation_duration:1/(2*Ti);
% plot(f,F);