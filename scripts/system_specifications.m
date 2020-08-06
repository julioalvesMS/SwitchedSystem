%% System Specifications

% Simulation specifications
Ti=1e-7; % [s] - Simulation pace
Tlog = 1e-5;
pwm_period = 1/20e3; % [s] PWM period
pwm_sample_time = Ti;
Tref = 1e-3;

%current_correction_start = 0.3;
current_correction_start = 0.2;

reference_start_time = 0.005;

% Circuit specifications
Ro = 96.8; % [Ohm] - Load Resistance
Co = 2250e-6; % [F] - Output Capacitance
Rc = 0.04;
R  = 0.49; % [Ohm] - Converter Resistance
L  = 1.981e-3; % [H] - Converter Indutance

tau = Ro*Co;

% Input voltage
Vs = 66; % [V]

% System starts discharged
x0 = [0; 0];

% Main control specifications
% Tsw = 5e-5; % [s] - Modern controller maximum switching period
Ts = 1/40e3; % [s] - Modern controller maximum switching period
%Ts = Ti;  % [s] - Modern controller maximum switching period
Fs = 1/Ts;   % [Hz] - Modern controller maximum switching frequency


% Parameters for dynamic reference
variable_reference_step_period = 0.5;


s = tf('s');
F = 1/(tau*s+1);
Fd = c2d(F, Tref, 'tustin');
[NFd, DFd] =  tfdata(Fd, 'v');