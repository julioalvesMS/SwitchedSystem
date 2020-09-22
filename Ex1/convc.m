function [ Ql ] = convc( Q, lmb )
%CONVC calculates the convex combination Ql of the polytopic Q cell of
%matrices with coefficients lmb

N = length(lmb);
Ql = Q{1}*0;
for i =1:N
    Ql= Q{i}*lmb(i)+Ql;
end

end

