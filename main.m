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
disturbance_Ro_enable = false;

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

% Desired DC-DC converter to use
% Options can be found in the system directory:
%   buck
%   boost
%   buck_boost
%   buck_boost_non_inverting
circuit = buck_boost(R, Ro, Co, L);


%test_voltages = circuit.test_voltages;
test_voltages = circuit.single_voltage;
%test_voltages = 300;

% test_voltages = [190];

simulation_duration = 1;


%% Measurements

opt_measurement_frequency = true;
opt_measurement_efficiency = false;
opt_measurement_clock = true;


mu = -1
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

%% Lambdas to simulate

if sys.N == 2
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
            model = 'sim_converter.slx';
        else
            model = 'sim_discrete.slx';
        end
end

% Number of simulations to run
Ns = size(lambdas, 1);

bar = waitbar(0, 'Preparing simulation', 'name', 'Simulating');

try
    load_system(model);
    
    run comment_simulink

    for i=Ns:-1:1
        
        Vref = test_voltages(i);

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
        tic
        sim(model, simulation_duration);
        toc

        % Store only samples of the data, this will be made in order to save
        % memory use
        sim_out(i).IL = get_logged_data(logsout, 'IL', plot_compression_rate);
        sim_out(i).Vout = get_logged_data(logsout, 'Vout', plot_compression_rate);
        sim_out(i).xe = get_logged_data(logsout, 'xe', plot_compression_rate);
        sim_out(i).Vref = get_logged_data(logsout, 'Vref', plot_compression_rate);
        sim_out(i).F = get_logged_data(logsout, 'F', plot_compression_rate);
        sim_out(i).Eff = get_logged_data(logsout, 'Eff', plot_compression_rate);
    end
    
    uncomment_blocks(model)
    
catch exception
    close(bar);
    rethrow(exception);
end
close(bar);

%%

if sim_out(1).F.Data
    dt = 1e-6;
    F = [];
    FV = [];

    for t=0.9:0.5:max(sim_out(i).F.Time)
        data = getsampleusingtime(sim_out(i).F, t-dt, t+dt);
        Rdata = getsampleusingtime(sim_out(i).Vref, t-dt, t+dt);

        F(end+1) = mean(data.Data);
        FV(end+1) = mean(Rdata.Data);
    end

    figure
    plot(FV, F, 'd-');
end
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


