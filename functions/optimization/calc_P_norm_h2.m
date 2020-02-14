function [P,K,W] = calc_P_norm_h2(dsys)

    A = dsys.A;
    B = dsys.B;
    E = dsys.C;
    F = dsys.D;
    
    [P,K,W] = solve_P_lmi(A, B, E, F);
end



function [P,K,Wo] = solve_P_lmi(A, B, E, F)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nx = size(A,1);
    nw = size(B,2);
    nu = size(B,2);
    %%%%%%%%%%%%%%%%%%%%%%%% NORMA H2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    % Descreve a LMI a ser projetada
    setlmis([])

    % declaracao de variveis
    %===============
    S = lmivar(1,[nx 1]);
    Y = lmivar(2,[nu nx]);
    W = lmivar(1,[nw 1]);
    %===============

    ct = -1;
    lmiterm([ct,1,1,S],1,1);
    lmiterm([ct,2,1,S],A,1);
    lmiterm([ct,2,1,Y],-B,1);
    lmiterm([ct,2,2,S],1,1);
    lmiterm([ct,3,1,S],E,1);
    lmiterm([ct,3,1,Y],-F,1);
    lmiterm([ct,3,3,0],1);

    ct = ct-1;
    lmiterm([ct,1,1,W],1,1);
    lmiterm([ct,2,1,0],B);
    lmiterm([ct,2,2,S],1,1);
    lmiterm([ct,3,3,0],1);

    lmisys = getlmis;


    % Declaracao funcao objetivo
    options = [1e-9,2000,0,200,1];
    %===============

    np = decnbr(lmisys);

    c = zeros(np,1);
    for i=1:np
        Wi = defcx(lmisys,i,W);
        c(i) = trace(Wi);
    end

    [copt,xopt] = mincx(lmisys,c,options);


    if (isempty(copt))
        ME = MException('LmiException:noPossibleResult', ...
        'A LMI não encontrou nenhuma resposta para P');
        throw(ME)
    end

    So = dec2mat(lmisys,xopt,S);
    Yo = dec2mat(lmisys,xopt,Y);
    Wo = dec2mat(lmisys,xopt,W);
    
    P = inv(So);
    K = Yo/So;
end