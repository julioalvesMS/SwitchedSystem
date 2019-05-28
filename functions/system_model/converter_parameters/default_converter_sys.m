function sys = default_converter_sys(circuit)
    %% System Specifications

    % Circuit specifications
    R  = 2; % [Ohm] - Converter Resistance
    Ro = 50; % [Ohm] - Load Resistance
    Co = 470e-6; % [F] - Load Capacitance
    L  = 500e-6; % [H]

    % Input voltage
    Vs = 100; % [V]

    % System starts discharged
    x0 = [0; 0];

    %% Prepare Data

    % Get the converter Space State representation and store it in the sys
    % struct. The struct contains the fields
    %   A - Cell Array with each subsystem A matrix - sys.A{i}
    %   B - Cell Array with each subsystem B matrix - sys.B{i}
    %   C - Cell Array with each subsystem C matrix - sys.C{i}
    %   D - Cell Array with each subsystem D matrix - sys.D{i}
    %   Q - Cell Array with each subsystem Q matrix - sys.Q{i}
    %   N - Number of avaliable subsystems - sys.N
    sys = circuit.get_sys(R, Ro, Co, L);

    % Complete the sys struct with this data
    sys.U = Vs;
    sys.x0 = x0;
end

