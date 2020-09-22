function [Po, dsys,Wo, ko] = calc_sys_discrete_theorem_1_test_Wk(sys, dsys, lambdas)
%CALC_SYS_CHAINED_1 Calculate de desired mean system using lambda
%   Calculate de desigred mean system from the switched system
%   Returns the system, P matrix and the equilibrium points
    
    I = eye(length(dsys.A{1}));
    
    for i=1:dsys.N
        dsys.L{i} = [(dsys.A{i}-I)  dsys.B{i}];
    end
    
    volu = 0;

    for k=1:size(lambdas, 1)
        [P,W] = frank_wolfe(sys, dsys, lambdas, k);
        
        [Al, Bl] = calc_sys_lambda(dsys, lambdas(k,:));
        xe = -(Al-I)\Bl*dsys.U;
        [~,cr] = calc_sys_discrete_h(dsys, lambdas(k,:), P, xe, dsys.U);
        n = 2;
        vol{k} = (pi^(n/2)/(gamma((n/2)+1)))/sqrt(det(W{k}/cr));
        
        if vol{k} > volu
            volu = vol{k}
            Wo = W;
            Po = P;
            ko = k;
        end
            
        
    end
end

function [P, W] = frank_wolfe(sys, dsys, lambdas, k)
    
    E = 1e-3;
    I = eye(length(dsys.A{1}));
    
    [P,W] = solve_R_lmi(sys, dsys, lambdas, I, k);
        
    
    for i=1:100
        [R,S] = solve_R_lmi(sys, dsys, lambdas, W{k}, k);

        if abs(trace(W{k}\(S{k}-W{k}))) < E
            break;
        end

        fnc = @(alpha) (log(det(alpha*S{k}+(1-alpha)*W{k})));

        [~,alpha] = BuscaDicotomica(fnc, 0, 1, 1e20);

        P = alpha*R+(1-alpha)*P;
        
        for j=1:size(lambdas, 1)
            W{j} = alpha*S{j}+(1-alpha)*W{j};
        end
    end
end

function [Rs,Wk] = solve_R_lmi(sys, dsys, lambdas, Ws, j)

    I = eye(length(dsys.A{1}));

    nx = size(dsys.A{1}, 1);

    % Descreve a LMI a ser projetada
    setlmis([]);

    % declaracao de variveis
    P = lmivar(1,[nx 1]);
    
    for k=1:size(lambdas, 1)
        W{k} = lmivar(1,[nx 1]);
    end
    
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
        lmiterm([-ct,1,1,W{k}],-1,1)

%         % Ql > W
%         ct = newlmi;
%         lmiterm([ct,1,1,W],1,1)
%         lmiterm([-ct,1,1,0],Qlamb)

        % (ye ? xe)W(ye ? xe) < \sum ?i li P li,
        ct = newlmi;
        for i=1:length(lambda)
            lmiterm([-ct,1,1,P],lambda(i)*l{i}',l{i})
        end
        lmiterm([ct,1,1,W{k}],(ye-xe)',(ye-xe))
    
        % W > 0
        ct = newlmi;
        lmiterm([-ct,1,1,W{k}],1,1)
    end
    
    % P > 0
    ct = newlmi;
    lmiterm([-ct,1,1,P],1,1)
    
    lmisys = getlmis;


    % Declaracao funcao objetivo
    options = [1e-4,2000,1e9,200,1];
    %===============

    np = decnbr(lmisys);

    
    c = zeros(np,1);
    for i=1:np
        Wi = defcx(lmisys,i,W{j});
        c(i) = -trace(Ws\Wi);
    end

    [copt,xopt] = mincx(lmisys,c,options);

    if (isempty(copt))
        ME = MException('LmiException:noPossibleResult', ...
        'A LMI não encontrou nenhuma resposta para P');
        throw(ME)
    end

    Rs = dec2mat(lmisys,xopt,P);
    
    for k=1:size(lambdas, 1)
        Wk{k} = dec2mat(lmisys,xopt,W{k});
    end
end