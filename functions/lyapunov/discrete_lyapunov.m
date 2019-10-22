function v = discrete_lyapunov(P, h, d, Xi)
    %DISCRETE_LYAPUNOV Summary of this function goes here
    %   Detailed explanation goes here
    
    v = [1 Xi'] * [d h'; h P] * [1; Xi];
end

