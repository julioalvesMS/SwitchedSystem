function sys = gss(A, B, C, D, Q, E, G, H)
%GSS General Space State system
    sys.A = A;
    sys.B = B;
    sys.C = C;
    sys.D = D;
    sys.Q = Q;
    sys.E = E;
    sys.G = G;
    sys.H = H;
    
    sys.N = length(A);
    sys.discrete = false;
end