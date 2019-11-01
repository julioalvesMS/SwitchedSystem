function [str] = gen_fnc_getD(class, dd1)
    y = {
"    double Buck::GetD(double P[SYSTEM_ORDER][SYSTEM_ORDER], double h[SYSTEM_ORDER])"
"    {"
"        double d = 0;"
""
"        switch(controlStrategy)"
"        {"
"        case CS_DISCRETE_THEOREM_1:"
"            //"
"            // " + class + " Converter - Discrete Rule 1"
"            //"
"            d = " + sprintf('%g;',dd1);
"            break;"
""
"        default:"
"            break;"
"        }"
""
"        return d;"
"    }"
    };
    str=sprintf('%s\n',y{:});
end

