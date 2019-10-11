function dsys = discrete_gss(sys, T)
    %DISCRETE_GSS Summary of this function goes here
    %   Detailed explanation goes here

    n = length(sys.A{1});
    
    dsys = sys;
    for i=1:length(sys.A)
        dsys.A{i} = expm(sys.A{i}*T);
        dsys.B{i} = (expm(sys.A{i}*T) - eye(n))/sys.A{i}*sys.B{i};
    end
        
end

