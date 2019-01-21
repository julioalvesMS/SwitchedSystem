function sys = chained_ss(sys, lambda)
%CHAINED_SS Summary of this function goes here
    
    nA = size(sys(1).A);
    nB = size(sys(1).B);
    nC = size(sys(1).C);
    nD = size(sys(1).D);

    A = zeros(nA);
    B = zeros(nB);
    C = zeros(nC);
    D = zeros(nD);
    
    for i=1:length(sys)
        A = A + lambda(i)*sys(i).A;
        B = B + lambda(i)*sys(i).B;
        C = C + lambda(i)*sys(i).C;
        D = D + lambda(i)*sys(i).D;
    end
    
    sys = gss(A, B, C, D);
end

