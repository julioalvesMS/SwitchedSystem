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

% Desired DC-DC converter to use
% Options can be found in the system directory:
%   buck
%   boost
%   buck_boost
circuit = boost(R, Ro, Co, L);

% Lambda used to create the mean system with wich the controle will be
% designed
% Value must be a vector where:
%   sum(lambda) = 1
lambda = [0.5 0.5];

%% Prepare Data

run load_circuit_sys

%% System Equilibrium Points

[A, B] = calc_sys_lambda(sys, lambda);
C = sys.C{1};
D = sys.D{1};

xe = -A\B*Vs;
Ie = xe(1);
Ve = xe(2);

sys_ss = ss(A, B, C, D);

s = tf('s');
sys_tf = tf(sys_ss);

% buck
Kp = 0.123;
Ki = 28.9;
Kd = 0;
% Kp = 0.0062;
% Ki = 0.909;
% Kd = 0;

% boost
iKp = 0.0203;
iKi = 4.77;
iKd = 0;
vKp = 0.316;
vKi = 3.23;
vKd = 0;

N = 100;
Ci_pid = iKp + iKi/s + iKd * N/(1+N/s);

% Boost - Corrente/Duty
d = 0.5;
Gi = ((Ve/L)*s + (1-d)*Ie/(L*Co))/(s^2 - R*s/L + (1-d)^2/L*Co); 

% sisotool(Gi, Ci_pid);

Fi = feedback(Gi*Ci_pid,1)

Cv_pid = vKp + vKi/s + vKd * N/(1+N/s);

% sisotool(Fi, Cv_pid);


% Buck
Kb = Vs/(L*Co);
a1 = Ro/Co;
a2 = 1/(L*Co);
G = tf([0 Kb],[1 a1 a2])*Vs;

C_buck = Kp + Ki/s + Kd * N/(1+N/s);

sisotool(G, C_buck)




Cd_pid = c2d(Ci_pid, pwm_period);


[num, den] = tfdata(Cd_pid);
num = num{1};
den = den{1};


