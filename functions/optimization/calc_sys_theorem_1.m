function [P, xe] = calc_sys_theorem_1(sys, lambda)
%CALC_SYS_CHAINED_1 Calculate de desired mean system using lambda
%   Calculate de desigred mean system from the switched system
%   Returns the system, P matrix and the equilibrium points
    
    [Al, Bl, Ql] = calc_sys_lambda(sys, lambda);
    xe = -Al\Bl*sys.U;
    
    %P = lyap(Al', Ql);
    P = solve_P_lmi(sys, Al, Ql, xe);
end


function Po = solve_P_lmi(sys, Alamb, Qlamb, xe)

    x0 = sys.x0;
    
    nx = size(Alamb,1);

    % Descreve a LMI a ser projetada
    setlmis([]);

    % declaracao de variveis
    P = lmivar(1,[nx 1]);
    
    ct = newlmi;
    lmiterm([ct,1,1,P],Alamb',1,'s')
    lmiterm([ct,1,1,0],Qlamb)

    
    ct = newlmi;
    lmiterm([-ct,1,1,P],1,1)
    
    lmisys = getlmis;


    % Declaracao funcao objetivo
    options = [1e-4,2000,1e9,200,1];
    %===============

    np = decnbr(lmisys);

    c = zeros(np,1);
    for i=1:np
        Pi = defcx(lmisys,i,P);
        c(i) = (x0-xe)'*Pi*(x0-xe);
    end

    [copt,xopt] = mincx(lmisys,c,options);


    if (isempty(copt))
        ME = MException('LmiException:noPossibleResult', ...
        'A LMI não encontrou nenhuma resposta para P');
        throw(ME)
    end

    Po = dec2mat(lmisys,xopt,P);
end