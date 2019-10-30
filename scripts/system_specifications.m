%% System Specifications

% Simulation specifications
Ti=1e-6; % [s] - Simulation pace
pwm_period = 1e-4; % [s] PWM period
pwm_sample_time = 1e-6;

sensor_sample = 3e-5;

% Circuit specifications
% R  = 0.0105; % [Ohm] - Converter Resistance
% L  = 4.8e-6; % [H] - Converter Indutance
% Ro = 7.5; % [Ohm] - Load Resistance
% Co = 726e-6; % [F] - Output Capacitance
% Rc = 0.00015; % [Ohm] - Capcitor Resistance
R  = 0.4; % [Ohm] - Converter Resistance
L  = 1.954e-3; % [H] - Converter Indutance
Ro = 96.8; % [Ohm] - Load Resistance
Co = 2250e-6; % [F] - Output Capacitance
Rc = 0.00015; % [Ohm] - Capcitor Resistance

tau = Ro*Co;

% Input voltage
Vs = 65; % [V]

% System starts discharged
x0 = [0; 0];

% Main control specifications
% Tsw = 5e-5; % [s] - Modern controller maximum switching period
Tsw = -1; % [s] - Modern controller maximum switching period
Fsw = 1/Tsw; % [Hz] - Modern controller maximum switching frequency

Ts = 5e-5;

% Parameters for dynamic reference
Vo = 45;
ts = 0.2;

controller_enable_time = 0.2
rederence_pid_initial_condition = 65