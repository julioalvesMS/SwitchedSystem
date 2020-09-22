function [P, W] = calc_sys_discrete_theorem_1_range(sys, dsys, lambdas)
%CALC_SYS_CHAINED_1 Calculate de desired mean system using lambda
%   Calculate de desigred mean system from the switched system
%   Returns the system, P matrix and the equilibrium points
    
    fnc = @(beta) (calc_volume(sys, dsys, lambdas, beta));

    [~,beta] = BuscaDicotomica(fnc, 1, 1e5, 1e7);

    [P,W] = frank_wolfe(sys, dsys, lambdas, beta);
    
end

function V = calc_volume(sys, dsys, lambdas, beta)
    try
        [P,~] = frank_wolfe(sys, dsys, lambdas, beta);
    %     V = -trace(W\P);
        V = -log(det(P));
        if V < 0
            V = inf;
        end
    catch
        V = inf;
    end
end

function [P, W] = frank_wolfe(sys, dsys, lambdas, beta)
    
    E = 1e-3;
    I = eye(length(dsys.A{1}));
    
    [P,W] = solve_R_lmi(sys, dsys, lambdas, beta, I);
    
    for i=1:100
        [R,S] = solve_R_lmi(sys, dsys, lambdas, beta, P);

        if abs(trace(P\(R-P))) < E
            break;
        end

        fnc = @(alpha) (log(det(alpha*R+(1-alpha)*P)));

        [~,alpha] = BuscaDicotomica(fnc, 0, 1, 1e20);

        P = alpha*R+(1-alpha)*P;
        W = alpha*S+(1-alpha)*W;
    end
end

function [Rs,Ws] = solve_R_lmi(sys, dsys, lambdas, beta, Ps)

    I = eye(length(dsys.A{1}));

    nx = size(dsys.A{1}, 1);

    % Descreve a LMI a ser projetada
    setlmis([]);

    % declaracao de variveis
    P = lmivar(1,[nx 1]);
    W = lmivar(1,[nx 1]);
    
    for k=1:size(lambdas, 1)
        lambda = lambdas(k,:);
        [Alc, Blc] = calc_sys_lambda(sys, lambda);
        ye = -Alc\Blc*sys.U;

        [Al, Bl] = calc_sys_lambda(dsys, lambda);
        xe = -(Al-I)\Bl*dsys.U;
    
        for i=dsys.N:-1:1
            l{i} = (dsys.A{i}-I)*xe + dsys.B{i}*dsys.U;
        end

        % \sum ?i Ai P Ai ?P < ?W
        ct = newlmi;
        for i=1:length(lambda)
            lmiterm([ct,1,1,P],lambda(i)*dsys.A{i}',dsys.A{i})
        end
        lmiterm([ct,1,1,P],-1,1)
        lmiterm([-ct,1,1,W],-1,1)

        % \sum ?i Ai P Ai ?P < ?W
        ct = newlmi;
        lmiterm([-ct,1,1,0],1)
        for i=1:length(lambda)
            lmiterm([-ct,1,1,P],-lambda(i)*beta*l{i}',l{i})
            lmiterm([-ct,1,2,P],-lambda(i)*l{i}',dsys.A{i}/(eye(nx) - Al))
        end
        lmiterm([-ct,2,2,P],1,1)
        lmiterm([-ct,2,3,P],1,1)
        lmiterm([-ct,3,3,W],beta,1)

%         % Ql > W
%         ct = newlmi;
%         lmiterm([ct,1,1,W],1,1)
%         lmiterm([-ct,1,1,0],Qlamb)

        % (ye ? xe)W(ye ? xe) < \sum ?i li P li,
        ct = newlmi;
        for i=1:length(lambda)
            lmiterm([-ct,1,1,P],lambda(i)*l{i}',l{i})
        end
        lmiterm([ct,1,1,W],(ye-xe)',(ye-xe))
    end
    
    % P > 0
    ct = newlmi;
    lmiterm([-ct,1,1,P],1,1)
    
    % W > 0
    ct = newlmi;
    lmiterm([-ct,1,1,W],1,1)
    
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
    Ws = dec2mat(lmisys,xopt,W);
end