function [P, gamma] = calc_sys_discrete_theorem_1_test_1(sys, dsys, lambdas)
%CALC_SYS_CHAINED_1 Calculate de desired mean system using lambda
%   Calculate de desigred mean system from the switched system
%   Returns the system, P matrix and the equilibrium points

    [P,gamma] = solve_R_lmi(sys, dsys, lambdas);
end

function [Rs,gammas] = solve_R_lmi(sys, dsys, lambdas)

    I = eye(length(dsys.A{1}));
    
    E = [0.6  0
        0   1
        ];

    nx = size(dsys.A{1}, 1);

    % Descreve a LMI a ser projetada
    setlmis([]);

    % declaracao de variveis
    P = lmivar(1,[nx 1]);
    gamma = lmivar(1,[1 0]);
    
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
        lmiterm([ct,1,1,gamma],E'*E,1)
        
        ct = newlmi;
        for i=1:length(lambda)
            lmiterm([ct,1,1,P],lambda(i)*l{i}',l{i})
        end
        lmiterm([ct,1,1,0],-1)

%         % Ql > W
%         ct = newlmi;
%         lmiterm([ct,1,1,W],1,1)
%         lmiterm([-ct,1,1,0],Qlamb)

        % (ye ? xe)W(ye ? xe) < \sum ?i li P li,
        ct = newlmi;
        lmiterm([ct,1,1,gamma],(ye-xe)'*E',E*(ye-xe))
        lmiterm([ct,1,1,0],-1)
    end
    
    % P > 0
    ct = newlmi;
    lmiterm([-ct,1,1,P],1,1)
    
    % W > 0
    ct = newlmi;
    lmiterm([-ct,1,1,gamma],1,1)
    
    lmisys = getlmis;


    % Declaracao funcao objetivo
    options = [1e-4,2000,1e9,200,1];
    %===============

    np = decnbr(lmisys);

    c = zeros(np,1);
    for i=1:np
        gammai = defcx(lmisys,i,gamma);
        c(i) = -gammai;
    end

    [copt,xopt] = mincx(lmisys,c,options);

    if (isempty(copt))
        ME = MException('LmiException:noPossibleResult', ...
        'A LMI não encontrou nenhuma resposta para P');
        throw(ME)
    end

    Rs = dec2mat(lmisys,xopt,P);
    gammas = dec2mat(lmisys,xopt,gamma);
end