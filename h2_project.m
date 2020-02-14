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

load('data/CURVAS_02_04_compare.mat')
load('data/TESTE_02_06.mat')
load('data/TESTE_02_13.mat')

%% System Specifications

run system_specifications

%% Simulation Parametersr
opt_pwm = true;

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
opt_current_correction = false;

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

% Desired DC-DC converter to use
% Options can be found in the system directory:
%   buck
%   boost
%   buck_boost
circuit = buck(R, Ro, Co, L);

% Lambda used to create the mean system with wich the controle will be
% designed
% Value must be a vector where:
%   sum(lambda) = 1
lambda = [0.5 0.5];

%% Prepare Data

run load_circuit_sys

run circuit_disturbance

current_correction_gain =  circuit.current_correction_gain;


reference_pid_kp = circuit.reference_pid_kp;
reference_pid_ki = circuit.reference_pid_ki;

[reference_ve_limit_lower, reference_ve_limit_upper] = get_reference_ve_limits(circuit, Vs);

%% System Equilibrium Points

A = [-R/L  -1/L
      1/Co  -1/(Ro*Co)
];
B = [Vs/L; 0];
C = [1.8 1];
% C = [2 1];
%C = [sqrt(sys.Q{1}(1,1)) sqrt(sys.Q{1}(2,2))];
D = 0;

sys_ss = ss(A, B, C, D);

sys_tf = tf(sys_ss);


% Buck
Kb = Vs/(L*Co);
a1 = Ro/Co;
a2 = 1/(L*Co);
G = tf([0 Kb],[1 a1 a2])*Vs;
[num, den] = tfdata(sys_tf);
[As, Bs, Cs, Ds] = tf2ss(num{1}, den{1});
sys_h = ss(As, Bs, Cs, Ds);


[K,M] = project_state_feedback_h2(sys_ss, pwm_period)
% 
% syms Z taus
% [num, den] = tfdata(S*M);
% F_syms = poly2sym(cell2mat(num),Z)/poly2sym(cell2mat(den),Z);
% a = taus/1e16;
% b = taus/2e16;
% c = 2.5e-5/4;
% Mfs = (Z^2*(a+b)+Z*(-b-a*c))/(Z^2+Z*(-c-1)+c);
% Z = 1;
% 
% a = 1/1e5;
% b = -1/2;
% c = 1e-5;
% Mf = (1*z-1)/(z-1);
% [Mnum, Mden] = tfdata(Mf, 'v');
% eval(diff(Mfs*F_syms))
% 
% %step(Mf, 5)
% %step(Mf*((1.25e-05*z + 1.25e-05)/(z-1)), 5)
% 
% F = S*M;

Vref = 30;

lambdas = generate_lambda_voltage(sys, Vref);
[~, xe] = calc_sys_theorem_1(sys, lambdas(1,:));

Yref = C*xe;

%%
close all


sim('sim_feedback', 0.1)
i=1;
plot_compression_rate = 1
sim_out(i).IL = downsample_timeseries(logsout.get('IL').Values, plot_compression_rate);
sim_out(i).Vout = downsample_timeseries(logsout.get('Vout').Values, plot_compression_rate);
sim_out(i).xe = downsample_timeseries(logsout.get('xe').Values, plot_compression_rate);
sim_out(i).Vref = downsample_timeseries(logsout.get('Vref').Values, plot_compression_rate);

close all;

EXP_CONV = Buck_test;
EXP_CONV = Buck_H2;
plot_voltage_time(sim_out, circuit.name, image_folder);
hold on
plot(EXP_CONV.t*1e3 - 0.4, EXP_CONV.Vof - 0.37)

plot_current_time(sim_out, circuit.name, image_folder);
hold on
plot(EXP_CONV.t*1e3 - 0.4, EXP_CONV.IL + 0.3667)

%% Classic Control

Cp = circuit.pwm_pid_kp + circuit.pwm_pid_ki/s;

%sisotool(sys_ss, Cp)