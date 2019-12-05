function [y,x,k,data] = BuscaDicotomica(fnc, inferior, superior, precisao)
%BuscaDicotomica Summary of this function goes here
%   Detailed explanation goes here

    l = (superior-inferior)/precisao;
    %l = 1e-10;
    
    E = l/2 - l/1000;
    condicao_parada = 0;
    iteration_limit = 1e4;
    
    for k=0:iteration_limit
        
        x = (inferior+superior)/2;
        
        traceback.k(k+1) = k;
        traceback.x(k+1) = x;
        traceback.y(k+1) = fnc(x);
        
        
        if superior - inferior <= l
            condicao_parada = 1;
            break;
        end
        
        meio = (inferior+superior)/2;
        
        lambda = meio - E;
        mi = meio + E;
        
        fl = fnc(lambda);
        fm = fnc(mi);
        
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
    
	y = fnc(x);
    
    data.traceback = traceback;
    data.stop_condition = condicao_parada;
    data.iterations = k;
    data.iteration_limit = iteration_limit;
    data.method_name = 'Busca Dicotômica';
end

