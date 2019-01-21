function [A, B, C, D, Q] = gss2double(sys)
%GSS2DOUBLE General Space State system
    A = cat(3, sys.A);
    B = cat(3, sys.B);
    C = cat(3, sys.C);
    D = cat(3, sys.D);
    Q = cat(3, sys.Q);
end