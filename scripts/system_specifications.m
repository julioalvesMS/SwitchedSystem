%% System Specifications

% Simulation specifications
Ti=5e-7; % [s] - Simulation pace
pwm_period = 5e-6; % [Hz] PWM frequency
pwm_sample_time = 5e-7;

sensor_sample = 3e-5;

% Circuit specifications
R  = 0.0105; % [Ohm] - Converter Resistance
L  = 4.8e-6; % [H] - Converter Indutance
Ro = 7.5; % [Ohm] - Load Resistance
Co = 726e-6; % [F] - Output Capacitance
Rc = 0.00015; % [Ohm] - Capcitor Resistance

% Input voltage
Vs = 9; % [V]

% System starts discharged
x0 = [0; 0];

% Main control specifications
Fsw = 2.5e4; % [Hz] - Modern controller maximum switching frequency

% Parameters for dynamic reference
Vo = 270;
ts = 0.2;