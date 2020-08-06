function str = gen_fnc_getP(class, Pc1, Pc2, Pd1)
    y = {
"    void " + class + "::GetP(double P[SYSTEM_ORDER][SYSTEM_ORDER])"
"    {"
"        switch(controlStrategy)"
"        {"
"        case CS_CONTINUOUS_THEOREM_1:"
"            //"
"            // " + class + " Converter - Rule 1"
"            //"
"            P[0][0] = " + sprintf('%.10g;',Pc1(1,1));
"            P[0][1] = " + sprintf('%.10g;',Pc1(1,2));
"            P[1][0] = " + sprintf('%.10g;',Pc1(2,1));
"            P[1][1] = " + sprintf('%.10g;',Pc1(2,2));
"            break;"
""
"        case CS_CONTINUOUS_THEOREM_2:"
"            //"
"            // " + class + " Converter - Rule 2"
"            //"
"            P[0][0] = " + sprintf('%.10g;',Pc2(1,1));
"            P[0][1] = " + sprintf('%.10g;',Pc2(1,2));
"            P[1][0] = " + sprintf('%.10g;',Pc2(2,1));
"            P[1][1] = " + sprintf('%.10g;',Pc2(2,2));
"            break;"
""
"        case CS_DISCRETE_THEOREM_1:"
"            //"
"            // " + class + " Converter - Discrete Rule 1"
"            //"
"            P[0][0] = " + sprintf('%.10g;',Pd1(1,1)); 
"            P[0][1] = " + sprintf('%.10g;',Pd1(1,2)); 
"            P[1][0] = " + sprintf('%.10g;',Pd1(2,1)); 
"            P[1][1] = " + sprintf('%.10g;',Pd1(2,2)); 
"            break;"
""
"        default:"
"            break;"
"        }"
"    }"
    };

    str=sprintf('%s\n',y{:});
end

