function sys = gss(A, B, C, D, Q)
%GSS General Space State system
    sys.A = A;
    sys.B = B;
    sys.C = C;
    sys.D = D;
    sys.Q = Q;
    
    sys.N = length(A);
end