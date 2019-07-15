%% System Specifications

% Simulation specifications
Ti=1e-6; % [s] - Simulation pace
pwm_frequency = 1e5; % [Hz] PWM frequency
pwm_sample_frequency = 1e6;

% Circuit specifications
R  = 0.135; % [Ohm] - Converter Resistance
L  = 5e-3; % [H] - Converter Indutance
Ro = 96.8; % [Ohm] - Load Resistance
Co = 2250e-6; % [F] - Output Capacitance
Rc = 0.015; % [Ohm] - Capcitor Resistance

% Input voltage
Vs = 400; % [V]

% System starts discharged
x0 = [0; 0];

% Main control specifications
Fsw = 1e5; % [Hz] PWM frequency

% Parameters for dynamic reference
Vo = 270;
ts = 0.2;