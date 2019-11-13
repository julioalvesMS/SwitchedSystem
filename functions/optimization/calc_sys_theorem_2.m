function [P, xe] = calc_sys_theorem_2(sys, lambda)
%CALC_SYS_CHAINED_2 Calculate a P matrix to stable all subsystems
%   Returns the P matrix

% P Matrix from the article
%     P = [
%         L/2     0
%         0       Co/2
%     ];
    
    [Al, Bl, ~] = calc_sys_lambda(sys, lambda);
    xe = -Al\Bl*sys.U;
    
    P = solve_P_lmi(sys, xe);
end

function Po = solve_P_lmi(sys, xe)

    N = sys.N;
    A = sys.A;
    Q = sys.Q;
    x0 = sys.x0;
    
    nx = size(A{1},1);

    % Descreve a LMI a ser projetada
    setlmis([]);

    % declaracao de variveis
    P = lmivar(1,[nx 1]);
    
    for i=1:N
        ct = newlmi;
        lmiterm([ct,1,1,P],A{i}',1,'s')
        lmiterm([ct,1,1,0],Q{i})
        
    end
    
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