function [str] = gen_fnc_getStateFeedbackH2Controller(circuit, K, C, M)
    
    y = {
"    void " + circuit.class_name + "::GetStateFeedbackH2Controller(double K[2], double C[2], double* M)"
"    {"
"        K[0] = " + sprintf('%g;',K(1));
"        K[1] = " + sprintf('%g;',K(2));
""
"        C[0] = " + sprintf('%g;',C(1));
"        C[1] = " + sprintf('%g;',C(2));
""
"        (*M) = " + sprintf('%g;',M);
"    }"
    };
    str=sprintf('%s\n',y{:});

end

