%% Prepare Data

% Get the converter Space State representation and store it in the sys
% struct. The struct contains the fields
%   A - Cell Array with each subsystem A matrix - sys.A{i}
%   B - Cell Array with each subsystem B matrix - sys.B{i}
%   C - Cell Array with each subsystem C matrix - sys.C{i}
%   D - Cell Array with each subsystem D matrix - sys.D{i}
%   Q - Cell Array with each subsystem Q matrix - sys.Q{i}
%   N - Number of avaliable subsystems - sys.N
sys = circuit.get_sys();

% Complete the sys struct with this data
sys.U = Vs;
sys.x0 = x0;

dsys = discrete_gss(sys,Ts);

[pwm_limit_lower, pwm_limit_upper] = circuit.get_pwm_control_limits();

fnc_converter_Ie = circuit.get_converter_Ie_fnc();