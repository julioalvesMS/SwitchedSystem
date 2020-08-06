function P = calc_sys_theorem_2(sys)
%CALC_SYS_CHAINED_2 Calculate a P matrix to stable all subsystems
%   Returns the P matrix

% P Matrix from the article
%     P = [
%         L/2     0
%         0       Co/2
%     ];
    
    P = solve_P_lmi(sys);
end

function Po = solve_P_lmi(sys)

    N = sys.N;
    A = sys.A;
    Q = sys.Q;
    
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
        c(i) = trace(Pi);
    end

    [copt,xopt] = mincx(lmisys,c,options);

    if (isempty(copt))
        ME = MException('LmiException:noPossibleResult', ...
        'A LMI n�o encontrou nenhuma resposta para P');
        throw(ME)
    end

    Po = dec2mat(lmisys,xopt,P);
end