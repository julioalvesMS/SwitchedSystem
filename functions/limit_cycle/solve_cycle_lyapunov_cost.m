function [out] = solve_cycle_lyapunov_cost(A,mask,Q,x0,xe0)
    n = length(A{1});
    kappa = length(mask);
    setlmis([])

    P = cell(1, kappa);
    
    %Variables
    for i=1:kappa
        P{i} = lmivar(1,[n 1]);
        clmi = newlmi;
        lmiterm([-clmi,1,1,P{i}],1,1);
    end

    for i=1:kappa
        clmi = newlmi;
        %   Col1
        lmiterm([clmi,1,1,P{mod(i,kappa)+1}],A{mask(i)}',A{mask(i)});
        lmiterm([clmi,1,1,P{i}],-1,1);
        lmiterm([clmi,1,1,0],Q);
    end

    %Solution
    lmi_set = getlmis;

    n_dec = decnbr(lmi_set);
    c = zeros(n_dec,1);
    for i=1:n_dec
        Rn =  defcx(lmi_set,i,P{1});  
        c(i) = trace((x0-xe0)'*Rn*(x0-xe0));
    end


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
