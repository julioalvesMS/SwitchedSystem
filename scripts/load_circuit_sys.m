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

if (isa(circuit, 'boost'))
    Vi = Ro*Vs/(Ro + R);
else
    Vi = 0;
end
    

if(exist('opt_pwm','var') && opt_pwm)
    sensor_sample = pwm_period;
else
    sensor_sample = Ts;
end


dsys = discrete_gss(sys,Ts);

[pwm_limit_lower, pwm_limit_upper] = circuit.get_pwm_control_limits();

fnc_converter_Ie = circuit.get_converter_Ie_fnc();

reference_controller_initial = circuit.get_reference_initial(Vs);

% VOltage Filter
Bv = [0.000027213807988318872, 2*0.0000272138079883188725,0.000027213807988318872];
Av = [1,-1.985190657896261,0.9852995131282144 ];


%identify converter

is_buck = isa(circuit, 'buck');
is_boost = isa(circuit, 'boost');
is_buck_boost = isa(circuit, 'buck_boost');
is_buck_boost_3 = isa(circuit, 'buck_boost_non_inverting');