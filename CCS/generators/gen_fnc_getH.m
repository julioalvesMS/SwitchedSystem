function [str] = gen_fnc_getH(class, hd1)
    y = {
"    void " + class + "::GetH(double h[SYSTEM_ORDER])"
"    {"
"        switch(controlStrategy)"
"        {"
"        case CS_DISCRETE_THEOREM_1:"
"            //"
"            // " + class + " Converter - Discrete Rule 1"
"            //"
"            h[0] = " + sprintf('%g;',hd1(1));
"            h[1] = " + sprintf('%g;',hd1(2));
"            break;"
""
"        default:"
"            break;"
"        }"
"    }"
    };
    
    str=sprintf('%s\n',y{:});
end

