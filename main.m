%% Initial Setup
% clear; clc; close all;

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
Simulink.fileGenControl('set', 'CodeGenFolder', cache_folder);

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
opt_discrete = true;

% Desired Theorem to use
% Theorems defines
%   1 - Fixed Equilibrium
%   2 - Valiable Equilibrium
opt_theorem = 1;

% Project controller for a single point or for a range
% Options:
%   0 - Single Point
%   1 - Range of operation
opt_range_design = true;

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

% Update the equilibrium point from the system
% Options
%   0 - Use given equilibrium
%   1 - Update equilibrium based on given reference voltage
opt_update_equilibrium = true;

% Use a PI to determine and update the equilibrium point
% Options
%   0 - Don't use the PI
%   1 - Update the voltage from the equilibrium point using a PI controller
opt_equilibrium_controller = false;

% Use a PI to correct the equilibrium current estimation
% Options
%   0 - Don't use the PI
%   1 - Correct the current from the equilibrium point using a PI controller
opt_current_correction = true;

% Run with knowledge of only the equilibrium voltage.
% Uses a filter to estimate the equilibrium current.
% Options
%   0 - Don't use the filter
%   1 - Use the filter to estimate the equilibrium current
opt_partial_information = false;

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
disturbance_Ro_enable = true;

% Choose between a constant output voltage or one with a different profile
% Options
%   0 - Update reference according to the profile in the simulink
%   1 - Use constante reference
opt_variable_load = false;

% Simulate the dead time observed in real life switches
% Options
%   0 - Consider ideal switches
%   1 - Consider dead time
opt_dead_time = true;

% Simulate the dead time observed in real life switches
% Options
%   0 - Consider ideal switches
%   1 - Consider dead time
opt_mode_hopping = true;

% Simulate the dead time observed in real life switches
% Options
%   0 - Consider ideal switches
%   1 - Consider dead time
opt_sensor_noises = false;

% Desired DC-DC converter to use
% Options can be found in the system directory:
%   buck
%   boost
%   buck_boost
%   buck_boost_non_inverting
circuit = buck_boost(R, Ro, Co, L);


%test_voltages = circuit.test_voltages;
test_voltages = circuit.single_voltage;
% test_voltages = 5;

% test_voltages = [190];

mu = -1;

simulation_duration = 3;


% Ts = 1/1000e3; % [s] - Modern controller maximum switching period
% Fs = 1/Ts;   % [Hz] - Modern controller maximum switching frequency


%% Measurements

opt_measurement_frequency = false;
opt_measurement_efficiency = false;
opt_measurement_ripple = false;
opt_measurement_clock = false;
opt_measurement_error = false;


%% Prepare Data

current_correction_gain =  circuit.current_correction_gain;

pwm_pid_kp = circuit.pwm_pid_kp;
pwm_pid_ki = circuit.pwm_pid_ki;
        
pwm_pid_vc_vp = circuit.pwm_pid_vc_vp;
pwm_pid_vc_vi = circuit.pwm_pid_vc_vi;
pwm_pid_vc_cp = circuit.pwm_pid_vc_cp;
pwm_pid_vc_ci = circuit.pwm_pid_vc_ci;

reference_pid_kp = circuit.reference_pid_kp;
reference_pid_ki = circuit.reference_pid_ki;

[reference_ve_limit_lower, reference_ve_limit_upper] = get_reference_ve_limits(circuit, Vs);

run load_circuit_sys

run circuit_disturbance

%% Simulate Converter 

% Get model to simulate
switch(opt_model)
    case 1
        if opt_discrete == 0
            model = 'general_system';
        else
            model = 'discrete_general_system';
        end
    case 2
        if opt_discrete == 0
            model = 'sim_converter';
        else
            model = 'sim_discrete';
        end
end

% Number of simulations to run
Ns = size(test_voltages, 1);

bar = waitbar(0, 'Preparing simulation', 'name', 'Simulating');

try
    load_system(model);
    
    run comment_simulink

%     for i=Ns:-1:1
    for i=4
        
        Vref = test_voltages(1);
        
        if opt_range_design == 0
            lambdas = generate_lambda_voltage(sys, circuit.limit_cycle_voltage);
        else
            range = circuit.operation_range_voltage_min:1:circuit.operation_range_voltage_max;
            lambdas = generate_lambda_voltage(sys, range);
        end

        xe = calculate_equilibrium_point(circuit, Vs, circuit.limit_cycle_voltage);
        
        % Calculate the P matrix, as the equilibrium point. Calculation will be
        % in accordance with the chosen control theorem
        if opt_discrete == 0
            switch(opt_theorem)
                case 1
                    P = calc_sys_theorem_1_range(sys, lambdas);
                case 2
                    P = calc_sys_theorem_2(sys);
            end
        else
            [P, dsys] = calc_sys_discrete_theorem_1_range(sys, dsys, lambdas);
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
            SystemDataBus = create_bus_DiscreteSystemDataBus(Ad, Bd, P, Qd, sys.N, Ld);
        end
            

        % Update wait bar
        waitbar((Ns-i)/Ns, bar, sprintf('Simulation %i of %d',Ns-i+1, Ns));

        % Run simulation
        tic
        set_param(model,'SignalLogging','on')
        set_param(model, 'SignalLoggingName', 'logsout')
        set_param(model,'SimulationMode','rapid')
        sim(model, simulation_duration);
        toc

        % Store only samples of the data, this will be made in order to save
        % memory use
        sim_out(i).xe = get_logged_data(logsout, 'xe', plot_compression_rate);
        sim_out(i).IL = get_logged_data(logsout, 'IL', plot_compression_rate);
        sim_out(i).Vout = get_logged_data(logsout, 'Vout', plot_compression_rate);
        sim_out(i).Vref = get_logged_data(logsout, 'Vref', plot_compression_rate);
        sim_out(i).Verr = get_logged_data(logsout, 'Verr', plot_compression_rate);
        sim_out(i).F = get_logged_data(logsout, 'F', plot_compression_rate);
        sim_out(i).Eff = get_logged_data(logsout, 'Eff', plot_compression_rate);
        sim_out(i).Ripple = get_logged_data(logsout, 'Voltage Ripple', plot_compression_rate);
    end
    
    uncomment_blocks(model)
    
catch exception
    close(bar);
    rethrow(exception);
end
close(bar);
    plot_disturbance_voltage_time(sim_out, disturbance_Ro_time, circuit.name, image_folder);

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

%% Fim do script

if sim_out(1).Eff.Data
    sim_out(1).Eff.Data(end)
end
return
%%
close all

EXP_CONV = Buck_PWM;
plot_voltage_time(sim_out, circuit.name, image_folder);
hold on
plot(EXP_CONV.t*1e3 - 0.4, EXP_CONV.Vof - 0.38)

plot_current_time(sim_out, circuit.name, image_folder);
hold on
plot(EXP_CONV.t*1e3 - 0.4, EXP_CONV.IL + 0.3667)

%%
close all
load('matlab4_bb.mat')

figure
hold all
plot(th_2_dis_FV, th_2_dis_F/1e3, '+-black')
plot(th_2_dis_dd_FV, th_2_dis_dd_F/1e3, '+-r')
hold off
legend('Ideal switches', 'Switches with Dead-Time')
%set(legend, 'Position',[0.345238102475802 0.866269843540494 0.333928564190865 0.0869047596341088]);
ylabel('Switching frequency [kHz]')
xlabel('Output voltage [V]')
%set(gcf, 'Position',  [100, 100, 1600, 780])
set(gcf,'renderer','Painters')
saveas(gcf, strcat(image_folder, 'frequency_simulation_discrete_controller'), 'epsc');

figure
hold all
plot(th_2_con_FV, th_2_con_F/1e3, '+-black')
plot(th_2_con_dd_FV, th_2_con_dd_F/1e3, '+-r')
hold off
legend('Ideal switches', 'Switches with Dead-Time')
%set(legend, 'Position',[0.345238102475802 0.866269843540494 0.333928564190865 0.0869047596341088]);
ylabel('Switching frequency [kHz]')
xlabel('Output voltage [V]')
%set(gcf, 'Position',  [100, 100, 1600, 780])
set(gcf,'renderer','Painters')
saveas(gcf, strcat(image_folder, 'frequency_simulation_continuous_controller'), 'epsc');


