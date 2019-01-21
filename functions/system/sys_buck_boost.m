function sys = sys_buck_boost(R, R0, C0, L)
%SYS_BUCK_BOOST Space State from Buck-Boost DC-DC converter
    A1 = [
        -R/L  0
        0     -1/R0*C0
    ];


    A2 = [
        -R/L  -1/L
         1/C0 -1/R0*C0
    ];

    B1 = [1/L; 0];
    B2 = [0; 0];
    
    C1 = [1 0];
    C2 = C1;
    
    D1 = 0;
    D2 = D1;
    
    sys1 = gss(A1, B1, C1, D1);
    sys2 = gss(A2, B2, C2 ,D2);
    
    sys = [sys1 sys2];
end

