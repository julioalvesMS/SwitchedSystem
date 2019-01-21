clear; clc; close all;

%% Initial Setup
addpath(genpath('functions'))
addpath(genpath('models'))

%% System specifications

R  = 1; % [Ohm]
R0 = 1; % [Ohm]
C0 = 1; % [F]
L  = 1; % [H]

Vs = 1; % [V]

%% Common calculations
lambdas = generate_lambda_2d();
x0 = [1; 1];


%% Buck Converter 
sys = sys_buck(R, R0, C0, L);

N = length(sys);

lambdas = [0.5 0.5];
for i=1:size(lambdas, 1)
    sys_lambda = chained_ss(sys, lambdas(i,:));
    P = lyap(sys_lambda.A, sys_lambda.Q);
    xe = -sys_lambda.A\sys_lambda.B;
    
    [A, B, C, D, Q] = gss2double(sys);
    
    SystemDataBus = create_bus_SystemDataBus(A, B, P, Q, xe, N);

    sim('theorem_1.slx');
end


%% Boost Converter 
sys = sys_boost(R, R0, C0, L);

%% Buck-Boost Converter 
sys = sys_buck_boost(R, R0, C0, L);

