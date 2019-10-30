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

% for i=1:sys.N
%     n = length(sys.A{i});
%     
%     A = zeros(n+1);
%     A(1:n,1:n) = sys.A{i};
%     A(end,:) = [1/tau 0 -1/tau];
%     sys.A{i} = A;
%     
%     B = zeros(n+1, 1);
%     B(1:n) = sys.B{i};
%     sys.B{i} = B;
%     
%     C = zeros(1, n+1);
%     C(1:n) = sys.C{i};
%     sys.C{i} = C;
%     
%     Q = zeros(n+1);
%     Q(1:n,1:n) = sys.Q{i};
%     sys.Q{i} = Q;
% end
% 
% x0 = zeros(n+1,1);
% x0(1:n) = sys.x0;
% sys.x0 = x0;