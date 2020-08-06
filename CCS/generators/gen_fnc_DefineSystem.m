function [str] = gen_fnc_DefineSystem(sys)
    y = {
"    void DefineSystem()"
"    {"
"        SubSystem* subSys;"
    };
    
for i=1:sys.N
    id = sprintf("%d", i);
    aux = {
    ""
    "        //"
    "        // =============== Subsystem "+id+" ==============="
    "        //"
    "        subSys = &(system.subSystems[" + sprintf("%d", i-1) + "]);"
    ""
    "        //"
    "        // Subsystem "+id+" -- Matrix A"
    "        //"
    "        subSys->A[0][0] = " + sprintf('%.10g;',sys.A{i}(1,1));
    "        subSys->A[0][1] = " + sprintf('%.10g;',sys.A{i}(1,2));
    "        subSys->A[1][0] = " + sprintf('%.10g;',sys.A{i}(2,1));
    "        subSys->A[1][1] = " + sprintf('%.10g;',sys.A{i}(2,2));
    "        //"
    "        // Subsystem "+id+" -- Matrix B"
    "        //"
    "        subSys->B[0] = " + sprintf('%.10g;',sys.B{i}(1));
    "        subSys->B[1] = " + sprintf('%.10g;',sys.B{i}(2));
    "        //"
    "        // Subsystem "+id+" -- Matrix Q"
    "        //"
    "        subSys->Q[0][0] = " + sprintf('%.10g;',sys.Q{i}(1,1));
    "        subSys->Q[0][1] = " + sprintf('%.10g;',sys.Q{i}(1,2));
    "        subSys->Q[1][0] = " + sprintf('%.10g;',sys.Q{i}(2,1));
    "        subSys->Q[1][1] = " + sprintf('%.10g;',sys.Q{i}(2,2));
    ""
    };

    y = [y; aux];
end
    aux = {
"    }"
    };
    
    y = [y; aux];
    str=sprintf('%s\n',y{:});

end

