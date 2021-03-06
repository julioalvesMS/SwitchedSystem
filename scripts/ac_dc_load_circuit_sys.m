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
theta = 0;
sys.U = Vs*[cos(theta) cos(theta+2*pi/3) cos(theta+4*pi/3)]';
sys.x0 = x0;


if(exist('opt_pwm','var') && opt_pwm)
    sensor_sample = pwm_period;
else
    sensor_sample = Ts;
end


dsys = discrete_gss(sys,Ts);

[pwm_limit_lower, pwm_limit_upper] = circuit.get_pwm_control_limits();
