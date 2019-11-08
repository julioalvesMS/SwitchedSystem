function [str] = gen_fnc_getClassicController(circuit, T)
    
    s = tf('s');
    C = circuit.pwm_pid_kp + circuit.pwm_pid_ki/s;
    Cd = c2d(C, T);
    
    [num, den] = tfdata(Cd);
    num = num{1};
    den = den{1};
    
    y = {
"    void " + circuit.class_name + "::GetClassicController(double num[2], double den[2])"
"    {"
"        num[0] = " + sprintf('%g;',num(1));
"        num[1] = " + sprintf('%g;',num(2));
""
"        den[0] = " + sprintf('%g;',den(1));
"        den[1] = " + sprintf('%g;',den(2));
"    }"
    };
    str=sprintf('%s\n',y{:});

end

