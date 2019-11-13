function [str] = gen_fnc_DefineDiscreteSystem(dsys)
    y = {
"    void DefineDiscreteSystem()"
"    {"
"        SubSystem* subSys;"
    };
    
for i=1:dsys.N
    id = sprintf("%d", i);
    aux = {
    ""
    "        //"
    "        // =============== Subsystem "+id+" ==============="
    "        //"
    "        subSys = &(discreteSystem.subSystems[" + sprintf("%d", i-1) + "]);"
    ""
    "        //"
    "        // Subsystem "+id+" -- Matrix A"
    "        //"
    "        subSys->A[0][0] = " + sprintf('%g;',dsys.A{i}(1,1));
    "        subSys->A[0][1] = " + sprintf('%g;',dsys.A{i}(1,2));
    "        subSys->A[1][0] = " + sprintf('%g;',dsys.A{i}(2,1));
    "        subSys->A[1][1] = " + sprintf('%g;',dsys.A{i}(2,2));
    "        //"
    "        // Subsystem "+id+" -- Matrix B"
    "        //"
    "        subSys->L[0] = " + sprintf('%g;',dsys.L{i}(1));
    "        subSys->L[1] = " + sprintf('%g;',dsys.L{i}(2));
    "        //"
    "        // Subsystem "+id+" -- Matrix Q"
    "        //"
    "        subSys->Q[0][0] = " + sprintf('%g;',dsys.Q{i}(1,1));
    "        subSys->Q[0][1] = " + sprintf('%g;',dsys.Q{i}(1,2));
    "        subSys->Q[1][0] = " + sprintf('%g;',dsys.Q{i}(2,1));
    "        subSys->Q[1][1] = " + sprintf('%g;',dsys.Q{i}(2,2));
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

