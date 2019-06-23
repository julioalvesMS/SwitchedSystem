
%% Prepare Data

pwm_pid_kp = circuit.pwm_pid_kp;
pwm_pid_ki = circuit.pwm_pid_ki;
pwm_pid_kd = circuit.pwm_pid_kd;

reference_pid_kp = circuit.reference_pid_kp;
reference_pid_ki = circuit.reference_pid_ki;

[reference_ve_limit_lower, reference_ve_limit_upper] = get_reference_ve_limits(circuit, Vs);

run load_circuit_sys

run circuit_disturbance

%% Lambdas to simulate

lambdas = generate_lambda_voltage(sys, test_voltages);

%% Simulate Converter 

% Get model to simulate
switch(opt_model)
    case 1
        model = 'general_system.slx';
    case 2
        model = circuit.simulink;
end

% Number of simulations to run
Ns = size(lambdas, 1);

bar = waitbar(0, 'Preparing simulation', 'name', 'Simulating');

sim_param.SimulationMode = 'rapid';
sim_param.AbsTol         = '1e-5';

load_system(model);

sim_out = [];

for i=Ns:-1:1
    
    % Convert the truct used to represent the space state to double, so it
    % can be used in the simulink
    [A, B, C, D, Q] = gss2double(sys);
    
    % Calculate the P matrix, as the equilibrium point. Calculation will be
    % in accordance with the chosen control theorem
    switch(opt_theorem)
        case 1
            [P, xe] = calc_sys_theorem_1(sys, lambdas(i,:));
        case 2
            [P, xe] = calc_sys_theorem_2(sys, lambdas(i,:));
    end
    
    % Creates a bus, wich will be used in the simulink to simplify the
    % model
    SystemDataBus = create_bus_SystemDataBus(A, B, P, Q, sys.N);

    % Update wait bar
    waitbar((Ns-i)/Ns, bar, sprintf('Simulation %i of %d',Ns-i+1, Ns));
    
    % Run simulation
    sim(model, simulation_duration);
    
    % Store only samples of the data, this will be made in order to save
    % memory use
    sim_out(i).Iout = compress_data(logsout.get('Iout').Values, plot_compression_rate);
    sim_out(i).Vout = compress_data(logsout.get('Vout').Values, plot_compression_rate);
end
close(bar);

%% Analysis

plot_voltage_current(sim_out, circuit.name, image_folder);

plot_voltage_time(sim_out, circuit.name, image_folder);

plot_current_time(sim_out, circuit.name, image_folder);

if disturbance_Ro_enable || disturbance_Vin_enable
    plot_disturbance_voltage_time(sim_out, disturbance_Ro_time, circuit.name, image_folder);
end