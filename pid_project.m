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

Simulink.fileGenControl('set', 'CacheFolder', cache_folder);

%% Simulation Parametersr

% Desired DC-DC converter to use
% Options can be found in the system directory:
%   buck
%   boost
%   buck_boost
circuit = boost;

% Lambda used to create the mean system with wich the controle will be
% designed
% Value must be a vector where:
%   sum(lambda) = 1
lambda = [0.5 0.5];

sys = default_converter_sys(circuit);

%% System Equilibrium Points

[A, B] = calc_sys_lambda(sys, lambda);

C = sys.C{1};
D = sys.D{1};

sys_ss = ss(A, B, C, D);

s = tf('s');
sys_tf = tf(sys_ss);

Kp = 0.01;
Ki = 0.05 * Kp;
C_pid = Kp+Ki/s;

sisotool(sys_tf, C_pid);