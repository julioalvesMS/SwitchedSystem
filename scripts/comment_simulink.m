%% Controller
comment_blocks(model, 'PWM Control', opt_pwm)
comment_blocks(model, 'Switched Control', ~opt_pwm)

comment_blocks(model, 'Classic Voltage Current Controller', opt_pwm_current_controller)
comment_blocks(model, 'Classic Voltage Controller', ~opt_pwm_current_controller)

comment_blocks(model, 'Reference Update', opt_update_equilibrium)

comment_blocks(model, 'Equilibrium Control', opt_equilibrium_controller)

comment_blocks(model, 'Partial Information Filter', opt_partial_information)
comment_blocks(model, 'Equilibrium Current Determination', ~opt_partial_information)

comment_blocks(model, 'Current Correction', opt_current_correction)


%% Reference
comment_blocks(model, 'Constant Reference', opt_constant_reference)
comment_blocks(model, 'Variable Reference', ~opt_constant_reference)


%% Operarion
comment_blocks(model, 'Dead Time', opt_dead_time)

comment_blocks(model, 'Mode Hopping', opt_mode_hopping)

comment_blocks(model, 'Variable Load', opt_variable_load)
comment_blocks(model, 'Real Load', ~opt_variable_load)


%% Measurements
comment_blocks(model, 'Frequency Measurement', opt_measurement_frequency)
comment_blocks(model, 'Efficiency Measurement', opt_measurement_efficiency)
comment_blocks(model, 'Ripple Measurement', opt_measurement_ripple)
comment_blocks(model, 'Clock Measurement', opt_measurement_clock)
comment_blocks(model, 'Error Measurement', opt_measurement_error)
