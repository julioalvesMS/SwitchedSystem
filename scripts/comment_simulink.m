comment_blocks(model, 'PWM Control', opt_pwm)
comment_blocks(model, 'Switched Control', ~opt_pwm)

comment_blocks(model, 'Reference Update', opt_update_equilibrium)

comment_blocks(model, 'Equilibrium Control', opt_equilibrium_controller)

comment_blocks(model, 'Partial Information Filter', opt_partial_information)
comment_blocks(model, 'Equilibrium Current', ~opt_partial_information)

comment_blocks(model, 'Constant Reference', opt_constant_reference)
comment_blocks(model, 'Variable Reference', ~opt_constant_reference)