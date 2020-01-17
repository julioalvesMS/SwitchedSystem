function [out] = solve_cycle_lyapunov_hinf(A,mask,E,H,G)
    N = length(A);
    n = length(A{1});
    kappa = length(mask);
    setlmis([])

    P = cell(1, kappa);

    %Variables
    for i=1:kappa
        P{i} = lmivar(1,[n 1]);
    end
    rho= lmivar(1,[1 1]);

    for i=1:kappa
        clmi = newlmi;
        %   Col1
        lmiterm([-clmi,1,1,P{i}],1,1);
        lmiterm([-clmi,3,1,P{mod(i,kappa)+1}],1,A{mask(i)});
        lmiterm([-clmi,4,1,0],E{mask(i)});

        %   Col2
        lmiterm([-clmi,2,2,rho],1,1);
        lmiterm([-clmi,3,2,P{mod(i,kappa)+1}],1,H{mask(i)});
        lmiterm([-clmi,4,2,0],G{mask(i)});

        %   Col3
        lmiterm([-clmi,3,3,P{mod(i,kappa)+1}],1,1);

        %   Col4
        lmiterm([-clmi,4,4,0],1);
        for j= 1:N
            clmi = newlmi;
            lmiterm([clmi,1,1,P{mod(i,kappa)+1}],H{j}',H{j});
            lmiterm([clmi,1,1,0],G{j}'*G{j});
            lmiterm([clmi,1,1,rho],-1,1);

        end
    end

    %Solution
    lmi_set = getlmis;

    n_dec = decnbr(lmi_set);
    c = zeros(n_dec,1);
    c(decinfo(lmi_set,rho)) = 1;

    opts = [1e-4,2000,-1,10,1];
    [optval,xopt] = mincx(lmi_set,c,opts);
    if (isempty(xopt)) % No Feas ?
        out.feas = 0;
        return
    end
    
    out.feas = 1; % Feas

    out.optval = optval;
    %Variables
    for i =1:kappa
        out.P{i} = (dec2mat(lmi_set,xopt,P{i}));
    end
end