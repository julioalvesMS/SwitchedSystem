classdef simulation_config < matlab.mixin.Copyable
    
    properties
        opt_model
        opt_discrete
        opt_theorem
        opt_range_design
        opt_pwm
        opt_pwm_current_controller
        opt_update_equilibrium
        opt_equilibrium_controller
        opt_current_correction
        opt_partial_information
        opt_constant_reference
        variable_reference_step_period
        disturbance_Vin_enable
        disturbance_Ro_enable
        disturbance_Vin_time
        disturbance_Ro_time
        opt_variable_load
        opt_dead_time
        opt_mode_hopping
        opt_sensor_noises
        switching_period
        simulation_sample
        circuit
        test_voltages
        simulation_duration
        opt_measurement_frequency
        opt_measurement_efficiency
        opt_measurement_ripple
        opt_measurement_clock
        opt_measurement_error
        image_folder
        sim_out
    end
    
    methods
        function default_config = simulation_config()
            default_config.opt_model = 2;
            default_config.opt_discrete = false;
            default_config.opt_theorem = 2;
            default_config.opt_range_design = true;
            default_config.opt_pwm = false;
            default_config.opt_pwm_current_controller = false;
            default_config.opt_update_equilibrium = true;
            default_config.opt_equilibrium_controller = false;
            default_config.opt_current_correction = false;
            default_config.opt_partial_information = false;
            default_config.opt_constant_reference = true;
            default_config.disturbance_Vin_enable = false;
            default_config.disturbance_Ro_enable = false;
            default_config.disturbance_Vin_time = 1;
            default_config.disturbance_Ro_time = 1.5;
            default_config.opt_variable_load = false;
            default_config.opt_dead_time = true;
            default_config.opt_mode_hopping = true;
            default_config.switching_period = 1/40e3;
            default_config.simulation_sample = 1e-6;
            default_config.circuit = [];
            default_config.test_voltages = [];
            default_config.simulation_duration = 0.15;
            default_config.opt_measurement_frequency = false;
            default_config.opt_measurement_efficiency = false;
            default_config.opt_measurement_ripple = false;
            default_config.opt_measurement_clock = false;
            default_config.opt_measurement_error = false;
            default_config.image_folder = '';
            default_config.sim_out = [];
        end
    end
end

