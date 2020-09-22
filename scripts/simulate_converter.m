
%% Prepare Data

current_correction_gain =  circuit.current_correction_gain;

pwm_pid_kp = circuit.pwm_pid_kp;
pwm_pid_ki = circuit.pwm_pid_ki;
        
pwm_pid_vc_vp = circuit.pwm_pid_vc_vp;
pwm_pid_vc_vi = circuit.pwm_pid_vc_vi;
pwm_pid_vc_cp = circuit.pwm_pid_vc_cp;
pwm_pid_vc_ci = circuit.pwm_pid_vc_ci;

reference_pid_kp = circuit.reference_pid_kp;
reference_pid_ki = circuit.reference_pid_ki;

[reference_ve_limit_lower, reference_ve_limit_upper] = get_reference_ve_limits(circuit, Vs);

run load_circuit_sys

%% Simulate Converter 

% Get model to simulate
switch(opt_model)
    case 1
        if opt_discrete == 0
            model = 'general_system';
        else
            model = 'discrete_general_system';
        end
    case 2
        if opt_discrete == 0
            model = 'sim_converter';
        else
            model = 'sim_discrete';
        end
end

% Number of simulations to run
Ns = length(test_voltages);

bar = waitbar(0, 'Preparing simulation', 'name', 'Simulating');

clear sim_out

try
    load_system(model);
    
    run comment_simulink

    for i=Ns:-1:1
        
        Vref = test_voltages(i);
        
        if opt_range_design == 0
            lambdas = generate_lambda_voltage(sys, circuit.single_voltage);
        else
            range = circuit.operation_range_voltage_min:5:circuit.operation_range_voltage_max;
            lambdas = generate_lambda_voltage(sys, range);
        end

        xe = calculate_equilibrium_point(circuit, Vs, circuit.limit_cycle_voltage);
        
        % Calculate the P matrix, as the equilibrium point. Calculation will be
        % in accordance with the chosen control theorem
        if opt_discrete == 0
            switch(opt_theorem)
                case 1
                    P = calc_sys_theorem_1_range(sys, lambdas);
                case 2
                    P = calc_sys_theorem_2(sys);
            end
        else
            switch(opt_theorem)
                case 1
                    [P, W] = calc_sys_discrete_theorem_1_range(sys, dsys, lambdas);
                case 2
                    [P, W] = calc_sys_discrete_theorem_1_test_W(sys, dsys, lambdas);
                case 3
                    [P, gamma] = calc_sys_discrete_theorem_1_test_1(sys, dsys, lambdas);
            end
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
            SystemDataBus = create_bus_DiscreteSystemDataBus(Ad, Bd, P, Qd, sys.N, Ld);
        end
            

        % Update wait bar
        waitbar((Ns-i)/Ns, bar, sprintf('Simulation %i of %d',Ns-i+1, Ns));

        % Run simulation
        set_param(model,'SignalLogging','on')
        set_param(model,'SignalLoggingName','logsout')
        set_param(model,'SimulationMode','rapid')
        set_param(model,'AccelVerboseBuild','off')
        sim(model, simulation_duration);
        
        if ~exist('logsout','var')
            logsout = tmp_raccel_logsout;
        end

        % Store only samples of the data, this will be made in order to save
        % memory use
        sim_out(i).xe = get_logged_data(logsout, 'xe', plot_compression_rate);
        sim_out(i).IL = get_logged_data(logsout, 'IL', plot_compression_rate);
        sim_out(i).Vout = get_logged_data(logsout, 'Vout', plot_compression_rate);
        sim_out(i).Vref = get_logged_data(logsout, 'Vref', plot_compression_rate);
        sim_out(i).Verr = get_logged_data(logsout, 'Verr', plot_compression_rate);
        sim_out(i).F = get_logged_data(logsout, 'F', plot_compression_rate);
        sim_out(i).Eff = get_logged_data(logsout, 'Eff', plot_compression_rate);
        sim_out(i).Ripple = get_logged_data(logsout, 'Voltage Ripple', plot_compression_rate);
    end
    
    uncomment_blocks(model)
    
catch exception
    close(bar);
    rethrow(exception);
end
close(bar);

%% Analysis

% plot_voltage_current(sim_out, circuit.name, image_folder);
% 
% plot_voltage_time(sim_out, circuit.name, image_folder);
% 
% plot_current_time(sim_out, circuit.name, image_folder);
% 
% if disturbance_Ro_enable == 1
%     plot_disturbance_voltage_time(sim_out, disturbance_Ro_time, circuit.name, image_folder);
% end
% 
% if disturbance_Vin_enable == 1
%     plot_disturbance_voltage_time(sim_out, disturbance_Vin_time, circuit.name, image_folder);
% end
