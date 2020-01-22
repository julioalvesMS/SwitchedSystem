comment_blocks(model, 'PWM Control', opt_pwm)
comment_blocks(model, 'Switched Control', ~opt_pwm)

comment_blocks(model, 'Classic Voltage Current Controller', opt_pwm_current_controller)
comment_blocks(model, 'Classic Voltage Controller', ~opt_pwm_current_controller)

comment_blocks(model, 'Reference Update', opt_update_equilibrium)

comment_blocks(model, 'Equilibrium Control', opt_equilibrium_controller)

comment_blocks(model, 'Partial Information Filter', opt_partial_information)
comment_blocks(model, 'Equilibrium Current Determination', ~opt_partial_information)

comment_blocks(model, 'Current Correction', opt_current_correction)

comment_blocks(model, 'Constant Reference', opt_constant_reference)
comment_blocks(model, 'Variable Reference', ~opt_constant_reference)