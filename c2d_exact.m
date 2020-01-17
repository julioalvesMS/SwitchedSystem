function [ Ad, bd ] = c2d_exact( A, b, T )
%C2D_EXACT(A,B,T)
%   Outputs
%   - Ad, Bd

if(iscell(A))
    N = length(A); 
    n = size(A{1},1);
    m = size(b{1},2);
    for i=1:N
        M = expm([A{i} b{i};zeros(m,n) zeros(m)]*T);
        Ad{i} = M(1:n,1:n);
        bd{i} = M(1:n,n+1:n+m);
    end

else
    N=1;
    n = size(A,1);
    m = size(b,2);
    M = expm([A b;zeros(m,n) zeros(m)]*T);
    Ad = M(1:n,1:n);
    bd = M(1:n,n+1:n+m);
    

end


end


