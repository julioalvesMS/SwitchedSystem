function [str] = gen_fnc_getClassicVoltageCurrentController(circuit, T)
    
    s = tf('s');
    Cv = circuit.pwm_pid_vc_vp + circuit.pwm_pid_vc_vi/s;
    Cdv = c2d(Cv, T);
    Ci = circuit.pwm_pid_vc_cp + circuit.pwm_pid_vc_ci/s;
    Cdi = c2d(Ci, T);
    
    [vNum, vDen] = tfdata(Cdv);
    vNum = vNum{1};
    vDen = vDen{1};
    
    [iNum, iDen] = tfdata(Cdi);
    iNum = iNum{1};
    iDen = iDen{1};
    
    y = {
"    void " + circuit.class_name + "::GetClassicVoltageCurrentController(double vNum[2], double vDen[2], double iNum[2], double iDen[2])"
"    {"
"        vNum[0] = " + sprintf('%.10g;',vNum(1));
"        vNum[1] = " + sprintf('%.10g;',vNum(2));
""
"        vDen[0] = " + sprintf('%.10g;',vDen(1));
"        vDen[1] = " + sprintf('%.10g;',vDen(2));
""
"        iNum[0] = " + sprintf('%.10g;',iNum(1));
"        iNum[1] = " + sprintf('%.10g;',iNum(2));
""
"        iDen[0] = " + sprintf('%.10g;',iDen(1));
"        iDen[1] = " + sprintf('%.10g;',iDen(2));
"    }"
    };
    str=sprintf('%s\n',y{:});

end

