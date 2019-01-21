function sys = gss(A,B,C,D)
%GSS General Space State system
    sys.A = A;
    sys.B = B;
    sys.C = C;
    sys.D = D;
    sys.Q = C'*C;
end