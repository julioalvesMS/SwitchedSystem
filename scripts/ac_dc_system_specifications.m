%% System Specifications

% Simulation specifications
Ti=1e-6; % [s] - Simulation pace
pwm_period = 5e-5; % [s] PWM period
pwm_sample_time = 1e-6;
Tref = 1e-3;
Ts = -1;

current_correction_start = 0.2;

% Circuit specifications
Co = 2250e-6; % [F] - Output Capacitance
Rc = 0.04;
R  = 0.49; % [Ohm] - Converter Resistance
L  = 1.981e-3; % [H] - Converter Indutance

tau = R*Co;

% Input voltage
Vs = 220; % [V]

% System starts discharged
x0 = [100   0   0]';

% Main control specifications
% Tsw = 5e-5; % [s] - Modern controller maximum switching period
Tsw = -1; % [s] - Modern controller maximum switching period
Fsw = 1/Tsw; % [Hz] - Modern controller maximum switching frequency


% Parameters for dynamic reference
Vo = 45;
ts = 0.2;
