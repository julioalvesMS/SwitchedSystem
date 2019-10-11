function h = calc_sys_discrete_h(sys, lambda, P)
%calc_sys_discrete_h Summary of this function goes here

    sA = size(sys.A{1});
    sc = size(sys.B{1});

    Al = zeros(sA);
    cl = zeros(sc);
    
    I = eye(sA);
    
    for i=1:sys.N
        Al = Al + lambda(i)*sys.A{i};
        cl = cl + lambda(i)*sys.A{i}'*P*sys.B{i};
    end

    h = (I - Al')\cl;
end

