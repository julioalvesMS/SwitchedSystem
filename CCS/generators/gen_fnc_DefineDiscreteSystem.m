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
    "        subSys->A[0][0] = " + sprintf('%.10g;',dsys.A{i}(1,1));
    "        subSys->A[0][1] = " + sprintf('%.10g;',dsys.A{i}(1,2));
    "        subSys->A[1][0] = " + sprintf('%.10g;',dsys.A{i}(2,1));
    "        subSys->A[1][1] = " + sprintf('%.10g;',dsys.A{i}(2,2));
    "        //"
    "        // Subsystem "+id+" -- Matrix L"
    "        //"
    "        subSys->L[0][0] = " + sprintf('%.10g;',dsys.L{i}(1,1));
    "        subSys->L[0][1] = " + sprintf('%.10g;',dsys.L{i}(1,2));
    "        subSys->L[0][2] = " + sprintf('%.10g;',dsys.L{i}(1,3));
    "        subSys->L[1][0] = " + sprintf('%.10g;',dsys.L{i}(2,1));
    "        subSys->L[1][1] = " + sprintf('%.10g;',dsys.L{i}(2,2));
    "        subSys->L[1][2] = " + sprintf('%.10g;',dsys.L{i}(2,3));
    "        //"
    "        // Subsystem "+id+" -- Matrix Q"
    "        //"
    "        subSys->E[0][0] = " + sprintf('%.10g;',dsys.E{i}(1,1));
    "        subSys->E[0][1] = " + sprintf('%.10g;',dsys.E{i}(1,2));
    "        subSys->E[1][0] = " + sprintf('%.10g;',dsys.E{i}(2,1));
    "        subSys->E[1][1] = " + sprintf('%.10g;',dsys.E{i}(2,2));
    "        //"
    "        // Subsystem "+id+" -- Matrix H"
    "        //"
    "        subSys->H[0] = " + sprintf('%.10g;',dsys.H{i}(1));
    "        subSys->H[1] = " + sprintf('%.10g;',dsys.H{i}(2));
    "        //"
    "        // Subsystem "+id+" -- Matrix G"
    "        //"
    "        subSys->G[0] = " + sprintf('%.10g;',dsys.G{i}(1));
    "        subSys->G[1] = " + sprintf('%.10g;',dsys.G{i}(2));
    "        //"
    "        // Subsystem "+id+" -- Matrix Q"
    "        //"
    "        subSys->Q[0][0] = " + sprintf('%.10g;',dsys.Q{i}(1,1));
    "        subSys->Q[0][1] = " + sprintf('%.10g;',dsys.Q{i}(1,2));
    "        subSys->Q[1][0] = " + sprintf('%.10g;',dsys.Q{i}(2,1));
    "        subSys->Q[1][1] = " + sprintf('%.10g;',dsys.Q{i}(2,2));
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

