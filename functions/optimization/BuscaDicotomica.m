function [y,x,k,data] = BuscaDicotomica(fnc, inferior, superior, d, xi)
%BuscaDicotomica Summary of this function goes here
%   Detailed explanation goes here

    l = (superior-inferior)/1e10;
    %l = 1e-10;
    
    multidimensional = 0;
    if nargin >= 4
        multidimensional = 1;
    end
    
    E = l/2 - l/1000;
    condicao_parada = 0;
    iteration_limit = 1e4;
    
    for k=0:iteration_limit
        
        x = (inferior+superior)/2;
        
        traceback.k(k+1) = k;
        if multidimensional==1
            traceback.x(k+1,:) = (x*d+xi)';
            traceback.y(k+1) = fnc(x*d+xi);
        else
            traceback.x(k+1) = x;
            traceback.y(k+1) = fnc(x);
        end
        
        
        if superior - inferior <= l
            condicao_parada = 1;
            break;
        end
        
        meio = (inferior+superior)/2;
        
        lambda = meio - E;
        mi = meio + E;
        
        if multidimensional == 1
            fl = fnc(lambda*d + xi);
            fm = fnc(mi*d + xi);
        else
            fl = fnc(lambda);
            fm = fnc(mi);
        end
        
        if fl < fm
            superior = mi;
        elseif fl > fm
            inferior = lambda;
        else
            x = (mi+lambda)/2;
            superior = x;
            inferior = x;
        end
            
    end
    
    x = (inferior+superior)/2;
    
    if multidimensional == 1
        y = fnc(x*d + xi);
    else
        y = fnc(x);
    end
    
    data.traceback = traceback;
    data.stop_condition = condicao_parada;
    data.iterations = k;
    data.iteration_limit = iteration_limit;
    data.method_name = 'Busca Dicotômica';
end

