function [str] = gen_fnc_getCurrentCorrectionController(circuit, T)
    
    s = tf('s');
    CC = 1.5 + 100  /s;
    CCd = c2d(CC, T);
    
    DC = 1.5 + 100/s;
    DCd = c2d(DC, T);
    
    [Cnum, Cden] = tfdata(CCd, 'v');
    [Dnum, Dden] = tfdata(DCd, 'v');

    y = {
"    void " + circuit.class_name + "::GetCurrentCorrectionController(double num[2], double den[2])"
"    {"
"        switch(controlStrategy)"
"        {"
"        case CS_CONTINUOUS_THEOREM_1:"
"        case CS_CONTINUOUS_THEOREM_2:"
"            //"
"            // " + circuit.class_name + " Converter - Continuous"
"            //"
"            num[0] = " + sprintf('%.10g;',Cnum(1));
"            num[1] = " + sprintf('%.10g;',Cnum(2));
""
"            den[0] = " + sprintf('%.10g;',Cden(1));
"            den[1] = " + sprintf('%.10g;',Cden(2));
"            break;"
""
"        case CS_DISCRETE_THEOREM_1:"
"            //"
"            // " + circuit.class_name + " Converter - Discrete"
"            //"
"            num[0] = " + sprintf('%.10g;',Dnum(1));
"            num[1] = " + sprintf('%.10g;',Dnum(2));
""
"            den[0] = " + sprintf('%.10g;',Dden(1));
"            den[1] = " + sprintf('%.10g;',Dden(2));
"            break;"
""
"        default:"
"            break;"
"        }"
"    }"
    };

    str=sprintf('%s\n',y{:});

end

