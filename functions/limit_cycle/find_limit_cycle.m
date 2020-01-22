function [cycle] = find_limit_cycle(sys, kappa, cand, opt)

    A = sys.A;
    E = sys.E;
    H = sys.H;
    G = sys.G;
    b = sys.b;
    x0 = sys.x0;
    Q = sys.Q{1};
    N = length(sys.A);
    nw = size(sys.H{1},2);
    
    valopt= zeros(1,length(cand));
    out = cell(1,length(cand));
    
    switch opt
        case 1
            for i=1:length(cand)
                [out{i}]= solve_cycle_lyapunov_cost(A,cand{i}.mask,Q,x0,cand{i}.xe_h{1});
                if(out{i}.feas)
                    valopt(i) = out{i}.optval;
                else
                    valopt(i)=nan;
                end
            end
        case 2
            for i=1:length(cand)
                for j=1:N
                    ellm1{j} = A{j}*cand{i}.xe_h{kappa} -cand{i}.xe_h{1}+b{j};
                    L{j} = [ellm1{j*ones(nw,1)}];
                end
                [out{i}]= solve_cycle_lyapunov_h2(A,cand{i}.mask,E,H,G,L);
               % [out{i}]= opt_h2(A,cand{i}.mask,E,H,G,L,2);
                if(out{i}.feas)
                    valopt(i) = out{i}.optval;
                else
                    valopt(i)=nan;
                end
            end
        case 3
            for i=1:length(cand)
                [out{i}]= solve_cycle_lyapunov_hinf(A,cand{i}.mask,E,H,G);
                if(out{i}.feas)
                    valopt(i) = out{i}.optval;
                else
                    valopt(i)=nan;
                end
            end
    end
    if(all(isnan(valopt)))
        disp('infeas')
        return;
    end
    [~,idx] = min(valopt);
    
    xe_h = cand{idx}.xe_h;
    out = out{idx};
    for i=1:N
        for k=1:kappa
            ell{i,k} = A{i}*xe_h{k} -xe_h{mod(k,kappa)+1}+b{i};
        end
    end
    
    cycle = cand{idx};
    cycle.lyap = out;
    cycle.ell = ell;
    cycle.kappa = kappa;
    
end

