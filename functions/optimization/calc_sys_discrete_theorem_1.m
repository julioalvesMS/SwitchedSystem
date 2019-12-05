function [P, h, d, xe, dsys] = calc_sys_discrete_theorem_1(sys, dsys, lambda)
%CALC_SYS_CHAINED_1 Calculate de desired mean system using lambda
%   Calculate de desigred mean system from the switched system
%   Returns the system, P matrix and the equilibrium points
    
    E = 1e-3;
    I = eye(length(dsys.A{1}));
    
    [Alc, Blc] = calc_sys_lambda(sys, lambda);
    ye = -Alc\Blc*sys.U;
    
    [Al, Bl] = calc_sys_lambda(dsys, lambda);
    xe = -(Al-I)\Bl*dsys.U;
    
    for i=1:dsys.N
        dsys.L{i} = -(dsys.A{i}-I)*((Al-I)\Bl) + dsys.B{i};
        dsys.l{i} = (dsys.A{i}-I)*xe + dsys.B{i}*dsys.U;
    end
    
    P = I;
    for i=1:100
        R = solve_R_lmi(dsys, Al, lambda, P, xe, ye, 1e5);
        
        if abs(trace(P\(R-P))) < E
            break;
        end

        fnc = @(alpha) (log(det(alpha*R+(1-alpha)*P)));

        [~,alpha] = BuscaDicotomica(fnc, 0, 1);

        P = alpha*R+(1-alpha)*P;
    end
    h = calc_sys_discrete_h(dsys, lambda, P);
    d = h'/P*h;
end


function Rs = solve_R_lmi(sys, Alamb, lambda, Ps, xe, ye, beta)

    nx = size(Alamb,1);

    % Descreve a LMI a ser projetada
    setlmis([]);

    % declaracao de variveis
    P = lmivar(1,[nx 1]);
    W = lmivar(1,[nx 1]);
    
    % \sum ?i Ai P Ai ?P < ?W
    ct = newlmi;
    for i=1:length(lambda)
        lmiterm([ct,1,1,P],lambda(i)*sys.A{i}',sys.A{i})
    end
    lmiterm([ct,1,1,P],-1,1)
    lmiterm([-ct,1,1,W],-1,1)
    
    % \sum ?i Ai P Ai ?P < ?W
    ct = newlmi;
    lmiterm([-ct,1,1,0],1)
    for i=1:length(lambda)
        lmiterm([-ct,1,1,P],-lambda(i)*beta*sys.l{i}',sys.l{i})
        lmiterm([-ct,1,2,P],-lambda(i)*sys.l{i}',sys.A{i}/(eye(nx) - Alamb))
    end
    lmiterm([-ct,2,2,P],1,1)
    lmiterm([-ct,2,3,P],1,1)
    lmiterm([-ct,3,3,W],beta,1)

%     % Ql > W
%     ct = newlmi;
%     lmiterm([ct,1,1,W],1,1)
%     lmiterm([-ct,1,1,0],Qlamb)
    
    % P > 0
    ct = newlmi;
    lmiterm([-ct,1,1,P],1,1)
    % W > 0
    ct = newlmi;
    lmiterm([-ct,1,1,W],1,1)
    
    % (ye ? xe)W(ye ? xe) < \sum ?i li P li,
    ct = newlmi;
    for i=1:length(lambda)
        lmiterm([-ct,1,1,P],lambda(i)*sys.l{i}',sys.l{i})
    end
    lmiterm([ct,1,1,W],(ye-xe)',(ye-xe))
    
    lmisys = getlmis;


    % Declaracao funcao objetivo
    options = [1e-4,2000,1e9,200,1];
    %===============

    np = decnbr(lmisys);

    c = zeros(np,1);
    for i=1:np
        Pi = defcx(lmisys,i,P);
        c(i) = -trace(Ps\Pi);
    end

    [copt,xopt] = mincx(lmisys,c,options);

    if (isempty(copt))
        ME = MException('LmiException:noPossibleResult', ...
        'A LMI não encontrou nenhuma resposta para P');
        throw(ME)
    end

    Rs = dec2mat(lmisys,xopt,P);
end