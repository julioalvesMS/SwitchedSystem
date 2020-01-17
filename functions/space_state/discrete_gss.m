function dsys = discrete_gss(sys, T)
    %DISCRETE_GSS Summary of this function goes here
    %   Detailed explanation goes here

    n = length(sys.A{1});
    
    dsys = sys;
    for i=1:sys.N
        dsys.A{i} = expm(sys.A{i}*T);
        dsys.B{i} = (expm(sys.A{i}*T) - eye(n))/sys.A{i}*sys.B{i};
        dsys.b{i} = dsys.B{i}*dsys.U;
        dsys.L{i} = zeros(size(dsys.B{i}));
        dsys.l{i} = dsys.L{i};
    end
end

