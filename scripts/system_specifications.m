%% System Specifications

% Simulation specifications
Ti=1e-6; % [s] - Simulation pace
pwm_frequency = 1e5; % [H]
pwm_sample_frequency = 1e6;

% Circuit specifications
R  = 2; % [Ohm] - Converter Resistance
L  = 500e-6; % [H] - Converter Indutance
Ro = 50; % [Ohm] - Load Resistance
Co = 470e-6; % [F] - Output Capacitance

% Input voltage
Vs = 100; % [V]

% System starts discharged
x0 = [0; 0];