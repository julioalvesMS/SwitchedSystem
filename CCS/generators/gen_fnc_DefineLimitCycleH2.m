function [str] = gen_fnc_DefineLimitCycleH2(dsys, cycle)
    y = {
"    void DefineLimitCycleH2()"
"    {"
"        CycleStep* step;"
""
"        limitCycle.kappa = " + sprintf('%.10g;',cycle.kappa);
    };
    
for i=1:cycle.kappa
    id = sprintf("%d", i);
    aux = {
    ""
    "        //"
    "        // =============== Limit Cycle Step "+id+" ==============="
    "        //"
    "        step = &(limitCycle.cycleSteps[" + sprintf("%d", i-1) + "]);"
    ""
    "        //"
    "        // Cycle Step "+id+" -- Matrix P"
    "        //"
    "        step->P[0][0] = " + sprintf('%.10g;',cycle.lyap.P{i}(1,1));
    "        step->P[0][1] = " + sprintf('%.10g;',cycle.lyap.P{i}(1,2));
    "        step->P[1][0] = " + sprintf('%.10g;',cycle.lyap.P{i}(2,1));
    "        step->P[1][1] = " + sprintf('%.10g;',cycle.lyap.P{i}(2,2));
    "        //"
    "        // Cycle Step "+id+" -- Vector Xe"
    "        //"
    "        step->Xe[0] = " + sprintf('%.10g;',cycle.xe_h{i}(1));
    "        step->Xe[1] = " + sprintf('%.10g;',cycle.xe_h{i}(2));
    "        //"
    "        // Cycle Step "+id+" -- Matrix ell"
    "        //"
    };

    for j=1:dsys.N
        n = sprintf("%d", j-1);
        aux_ell = {
        "        step->ell["+n+"][0] = " + sprintf('%.10g;',cycle.ell{j,i}(1));
        "        step->ell["+n+"][1] = " + sprintf('%.10g;',cycle.ell{j,i}(2));
        };
        aux = [aux; aux_ell];
    end

    y = [y; aux];
end
    aux = {
    ""
    "    }"
    };
    
    y = [y; aux];
    str=sprintf('%s\n',y{:});

end

