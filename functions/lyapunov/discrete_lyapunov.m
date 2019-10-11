function v = discrete_lyapunov(h,P,csi)
    %DISCRETE_LYAPUNOV Summary of this function goes here
    %   Detailed explanation goes here
    
    d = 0;
    v = [1 csi'] * [d h'; h P] * [1; csi];
end

