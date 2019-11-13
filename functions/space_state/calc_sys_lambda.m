function [A, B, Q] = calc_sys_lambda(sys, lambda)
%SYS_CHAINED Summary of this function goes here
    
    sA = size(sys.A{1});
    sB = size(sys.B{1});
    sQ = size(sys.Q{1});

    A = zeros(sA);
    B = zeros(sB);
    Q = zeros(sQ);
    
    for i=1:sys.N
        A = A + lambda(i)*sys.A{i};
        B = B + lambda(i)*sys.B{i};
        Q = Q + lambda(i)*sys.Q{i};
    end
    
end

