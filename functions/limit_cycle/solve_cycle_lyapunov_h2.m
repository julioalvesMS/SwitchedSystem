function [out] = solve_cycle_lyapunov_h2(A,mask,E,H,G,L,mo)
    
    if(~exist('mo','var'))
        mo = 0;
    end
    
    N = length(A);
    n = length(A{1});
    nw = size(H{1},2);
    kappa = length(mask);
    setlmis([])

    P = cell(1, kappa);
    Q = cell(1, N);

    %Variables
    for i=1:kappa
        P{i} = lmivar(1,[n 1]);
    end
    for i=1:N
        Q{i} = lmivar(1,[nw 1]);
    end

    for i = 1:(kappa)
        clmi = newlmi;
        %   Col1
        lmiterm([clmi,1,1,P{mod(i,kappa)+1}],A{mask(i)}',A{mask(i)});
        lmiterm([clmi,1,1,P{i}],-1,1);
        lmiterm([clmi,1,1,0],E{mask(i)}'*E{mask(i)});
    end
    
    for i=1:N
        clmi = newlmi;
        lmiterm([-clmi,1,1,Q{i}],1,1);
        lmiterm([-clmi,2,1,P{1}],1,L{i}+H{i});
        lmiterm([-clmi,3,1,0],G{i});

        lmiterm([-clmi,2,2,P{1}],1,1);
        lmiterm([-clmi,3,3,0],1);
    end


    %Solution
    lmi_set = getlmis;

    n_dec = decnbr(lmi_set);

    xopt_ = cell(1, N);
    optval_ = inf(1, N);
    
    for m=1:N
        c = zeros(n_dec,1);
        for i=1:n_dec
            Rn =  defcx(lmi_set,i,Q{m});  
            c(i) = trace(Rn);
        end

        opts = [1e-4,2000,-1,10,1];
        [optval,xopt] = mincx(lmi_set,c,opts);
        if (isempty(xopt)) % No Feas ?
            out.feas = 0;
            return
        end
        xopt_{m} = xopt;
        optval_(m)=optval;
    end
    if(mo==0)
        [optval,m] = min(optval_);
    else
        m = mo;
        optval = optval_(mo);
    end

    xopt= xopt_{m};
    out.feas = 1; % Feas

    out.optval = optval;
    out.m=m;
    %Variables
    for i =1:kappa
        out.P{i} = (dec2mat(lmi_set,xopt,P{i}));
    end
end