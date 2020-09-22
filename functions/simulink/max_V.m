function [rs,betas] = max_V(W,P,h)
    % Descreve a LMI a ser projetada
    setlmis([])

    % Declaração de variáveis
    %===============
    r = lmivar(1,[1 0]);
    beta  = lmivar(1,[1 0]);
    %===============

    % Descrição das LMIs
    %===============
    ct = 0;

    ct = ct+1;
    lmiterm([-ct,1,1,r],1,1);
    lmiterm([-ct,1,1,beta],-1,1);
    lmiterm([-ct,2,1,0],h);
    lmiterm([-ct,2,2,0],P);
    lmiterm([-ct,3,2,0],P);
    lmiterm([-ct,3,3,beta],1,W);

    ct = ct+1;
    lmiterm([-ct,1,1,beta],1,1);

    lmisys = getlmis;
    %===============

    % Declaração função objetivo
    options = [1e-4,2000,1e9,200,0];
    %===============
    np = decnbr(lmisys);

    c = zeros(np,1);
    c(decinfo(lmisys,r)) = 1;

    [copt,xopt] = mincx(lmisys,c,options);

    if (isempty(copt))
        rs = NaN;
        betas = NaN;
    else
        rs = dec2mat(lmisys,xopt,r);
        betas = dec2mat(lmisys,xopt,beta);
    end
end
