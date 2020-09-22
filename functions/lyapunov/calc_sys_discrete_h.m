function [h, cr] = calc_sys_discrete_h(sys, lambda, P, xe, u)
%calc_sys_discrete_h Summary of this function goes here

    sA = size(sys.A{1});
    sc = size(sys.B{1});
    sr = size(sys.B{1});

    Al = zeros(sA);
    cl = zeros(sc);
    cr = zeros(1);
    
    I = eye(sA);
    
    for i=1:sys.N
        Al = Al + lambda(i)*sys.A{i};
        ell = sys.L{i}*[xe' u]';
        cl = cl + lambda(i)*sys.A{i}'*P*ell;
        cr = cr + lambda(i)*ell'*P*ell;
    end

    h = (I - Al')\cl;
end

