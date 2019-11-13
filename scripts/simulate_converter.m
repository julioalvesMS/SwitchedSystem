
%% Prepare Data

pwm_pid_kp = circuit.pwm_pid_kp;
pwm_pid_ki = circuit.pwm_pid_ki;
pwm_pid_kd = circuit.pwm_pid_kd;

reference_pid_kp = circuit.reference_pid_kp;
reference_pid_ki = circuit.reference_pid_ki;

[reference_ve_limit_lower, reference_ve_limit_upper] = get_reference_ve_limits(circuit, Vs);

run load_circuit_sys

%% Lambdas to simulate

lambdas = generate_lambda_voltage(sys, test_voltages);

%% Simulate Converter 

% Get model to simulate
switch(opt_model)
    case 1
        if opt_discrete == 0
            model = 'general_system.slx';
        else
            model = 'discrete_general_system.slx';
        end
    case 2
        if opt_discrete == 0
            model = circuit.simulink;
        else
            model = circuit.discrete_simulink;
        end
end

% Number of simulations to run
Ns = size(lambdas, 1);

bar = waitbar(0, 'Preparing simulation', 'name', 'Simulating');

clear sim_out

try
    load_system(model);

    for i=Ns:-1:1

        % Calculate the P matrix, as the equilibrium point. Calculation will be
        % in accordance with the chosen control theorem
        if opt_discrete == 0
            switch(opt_theorem)
                case 1
                    [P, xe] = calc_sys_theorem_1(sys, lambdas(i,:));
                case 2
                    [P, xe] = calc_sys_theorem_2(sys, lambdas(i,:));
            end
        else
            [P, h, d, xe, dsys] = calc_sys_discrete_theorem_1(sys, dsys, lambdas(i,:));
        end

        % Convert the truct used to represent the space state to double, so it
        % can be used in the simulink
        [A, B, C, D, Q] = gss2double(sys);          % Continuos
        [Ad, Bd, Cd, Dd, Qd, Ld] = gss2double(dsys);    % Discrete


        % Creates a bus, wich will be used in the simulink to simplify the
        % model
        if opt_discrete == 0
            SystemDataBus = create_bus_SystemDataBus(A, B, P, Q, sys.N);
        else
            SystemDataBus = create_bus_DiscreteSystemDataBus(Ad, Bd, P, Qd, sys.N, h, Ld);
        end
            

        % Update wait bar
        waitbar((Ns-i)/Ns, bar, sprintf('Simulation %i of %d',Ns-i+1, Ns));

        % Run simulation
        sim(model, simulation_duration);

        % Store only samples of the data, this will be made in order to save
        % memory use
        sim_out(i).IL = downsample(logsout.get('IL').Values, plot_compression_rate);
        sim_out(i).Vout = downsample(logsout.get('Vout').Values, plot_compression_rate);
        sim_out(i).xe = downsample(logsout.get('xe').Values, plot_compression_rate);
        sim_out(i).Vref = downsample(logsout.get('Vref').Values, plot_compression_rate);
    end
catch exception
    close(bar);
    rethrow(exception);
end
close(bar);

%% Analysis

plot_voltage_current(sim_out, circuit.name, image_folder);

plot_voltage_time(sim_out, circuit.name, image_folder);

plot_current_time(sim_out, circuit.name, image_folder);

if disturbance_Ro_enable == 1
    plot_disturbance_voltage_time(sim_out, disturbance_Ro_time, circuit.name, image_folder);
end

if disturbance_Vin_enable == 1
    plot_disturbance_voltage_time(sim_out, disturbance_Vin_time, circuit.name, image_folder);
end