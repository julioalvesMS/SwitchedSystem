function [C, kappa] = find_cycles(sys, xe, Gamma, K)
    C = {};
    
    if(~exist('K','var')|| isempty(K))
        K = 1:15;
    end
    
    for kappa=K
        C = find_cycles_kappa(sys, xe, Gamma, kappa);
        if ~isempty(C)
            break;
        end
    end
end

function C = find_cycles_kappa(sys, xe, Gamma, kappa)
    
    GAMMA = {Gamma};
    GX_ast = {(GAMMA{1}*xe)'};

    A = sys.A;
    b = sys.b;

    N = length(sys.A);
    n = length(sys.A{1});

    v_aux = cell(1, kappa);
    for i=1:kappa
        v_aux{i} = 1:N;
    end
    mask = combvec(v_aux{:});

    minusmat = -eye(n*kappa);
    minusmat = minusmat([(n+1):(n*kappa) 1:n],:);

    cand = {};
    for i = 1:size(mask,2)
        Atilde = blkdiag(A{mask(:,i)})+minusmat;
        btilde = cell2mat({b{mask(:,i)}}');
        xeb  = -Atilde\btilde;
        dists = reshape( blkdiag(GAMMA{ones(kappa,1)}) * xeb-[GX_ast{ones(kappa,1)}]', [length(GX_ast{1}), kappa]);

        val = (mean(sqrt(sum(dists.^2,1))));
        if val<1
            cand{end+1}.mask = mask(:,i);
            cand{end}.dist = val;
            cand{end}.singular = (min(abs(eig(Atilde)))<1e-10);
            for j=1:kappa
                cand{end}.xe_h{j} =  xeb((j-1)*n+(1:n));
            end
        end
    end
    C = cand;
end