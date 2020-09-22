function ex_find
ganhou =0;
try
    while(ganhou==0)
        F{1} = [-3 -6 3; 2 2 -3;  1.5582  0 -2];
        F{2} = [1 3 3 ;-0.1525 -3 -3; 0 0 -2];
        F{1} = [-3 -6 3; 2 2 -3;  1.6  0 -2];
        F{2} = [1 3 3 ;-0.2 -3 -3; 0 0 -2];
      %  F{1} = [-3 -6 3; 2 2 -3; rand*1.5+0.5 0 -2];
        %F{2} = [1 3 3 ;-rand -3 -3; 0 0 -2];
%         F{1} = [   2    -2     1;
%     -2    -2     2;
%     -2     0     1];
% F{2} = [-2    -2    -1;
%      1    -2    -1;
%     -2    -2    -1];
        g{1} = [0 1 0]';
        g{2} = [1 0 0]';
        N=2;
        n=3;
        T = 1%0.5;
        %ExLMdata
        [A,b] = c2d_exact(F,g,T);

        xe =  [-0.1550
    3.1035
   -6.7945]*0;
        for i =1:N
            ell{i} = (A{i}-eye(n))*xe+b{i};
        end

        span = 0:0.001:1;

        for i=1:length(span)
            lambda = [span(i) 1-span(i)];
            [outECC]= opt_ECC(A,ell,lambda);
            if(outECC.feas==1)
                outECC.lambda = lambda
                break;
            end

        end
        if(outECC.feas==1)
            fprintf('= ECC: %d \t LM: -\n',outECC.feas)
            continue;
        end

        span = 0:0.05:1;
        for i=1:length(span)
            for j=1:length(span)
                p = span(i);
                q = span(j);

                Mpi = [p  1-q;
                       1-p  q];
                [outLM]= opt_theorem1(A,ell,Mpi,[1 1]) ;
                if(outLM.feas==1)
                    outLM.Mpi = Mpi;
                    break;
                end
            end
            if(outLM.feas==1)
                break;
            end
        end
        if(outLM.feas==1)
            ganhou=1;
        end
        fprintf('= ECC: %d \t LM: %d\n',outECC.feas,outLM.feas)
    end
catch e
    keyboard
end
keyboard
end
%% opt_theorem1
function [out]= opt_theorem1(A,ell,Mpi,betas) 
N = length(A);
n = length(A{1});
setlmis([])

%Variables
W = lmivar(1,[n 1]);
h = lmivar(2,[n 1]);
for i =1:N
    P{i} = lmivar(1,[n 1]);
end
    
for i = 1:N
    %Sprocedure
    clmi = newlmi;
    %   Col1
    lmiterm([-clmi,1,1,P{i}],betas(i),1);
    for j=1:N
        lmiterm([-clmi,3,1,P{j}],betas(i)*Mpi(j,i),A{i});
    end
    lmiterm([-clmi,4,1,W],-1,1);
    %   Col2
    lmiterm([-clmi,2,2,0],1);
    for j=1:N
        lmiterm([-clmi,3,2,P{j}],betas(i)*Mpi(j,i),ell{i});
    end
    lmiterm([-clmi,4,2,h],-1,1);
    %   Col3
    for j=1:N
        lmiterm([-clmi,3,3,P{j}],Mpi(j,i),betas(i));
    end
    %   Col4
    lmiterm([-clmi,4,4,W],1,1);
end

%Solution
lmi_set = getlmis;
opts = [1e-4,2000,-1,10,1];
np = decnbr(lmi_set);

c = zeros(np,1);

[copt,xopt] = mincx(lmi_set,c,opts);

if (isempty(xopt)) % No Feas ?
    out.feas = 0;
    return
end

out.feas = 1; % Feas
out.W = dec2mat(lmi_set,xopt,W);
out.h = dec2mat(lmi_set,xopt,h);
for i=1:N
    out.P{i} = dec2mat(lmi_set,xopt,P{i});
end
%out.it = it;
end


%%  opt_IEEEtheorem1
function [out] = opt_IEEEtheorem1(A,ell,lambin)
n = size(A{1},1);
N = size(A,2);
setlmis([])

%Variables
P = lmivar(1,[n 1]);
W = lmivar(1,[n 1]);



ct = 1;
for i = 1:N
    lmiterm([ct,1,1,P],lambin(i)*A{i}',A{i});
end
lmiterm([ct,1,1,P],-1,1);
lmiterm([ct,1,1,W],1,1);

ct = ct+1;
for i=1:N
    lmiterm([ct,1,1,P],lambin(i)*ell{i}',ell{i});
end
lmiterm([ct,1,1,0],-1);

ct = ct+1;
lmiterm([-ct,1,1,P],1,1);

ct = ct+1;
lmiterm([-ct,1,1,W],1,1);

lmisys = getlmis;
%===============



options = [1e-4,2000,-1,200,1];
np = decnbr(lmisys);

c = zeros(np,1);

[copt,xopt] = mincx(lmisys,c,options);



if (isempty(copt))
    out.feas=0;
else
    out.feas=1;
    out.it= it;
    out.W = dec2mat(lmisys,xopt,W);
    out.P = dec2mat(lmisys,xopt,P);
    for i =1:N
        out.c{i} = A{i}'* out.P*ell{i};
        out.rho{i} = ell{i}'* out.P*ell{i};
    end
    out.h = (eye(n)-convc(A,lambin)')\convc(out.c,lambin);
end

end

%% opt_ECC
function [out]= opt_ECC(A,ell,lambda) 
N = length(A);
n = length(A{1});
setlmis([])

%Variaveis
W = lmivar(1,[n 1]);

P = lmivar(1,[n 1]);

clmi = newlmi;
for i = 1:N
    lmiterm([clmi,1,1,P],lambda(i)*A{i}',A{i});
end
lmiterm([clmi,1,1,P],-1,1);
lmiterm([clmi,1,1,W],1,1);

clmi = newlmi;
for i = 1:N
    lmiterm([-clmi,1,1,P],-lambda(i)*ell{i}',ell{i});
    lmiterm([-clmi,2,1,P],-lambda(i)*A{i}',ell{i});
end
    lmiterm([-clmi,1,1,0],1)
lmiterm([-clmi,2,2,W],1,1);

clmi = newlmi;
lmiterm([-clmi,1,1,P],1,1);


%Solucao
lmi_set = getlmis;

opts = [1e-4,2000,-1,10,1];

np = decnbr(lmi_set);

c = zeros(np,1);

[copt,xopt] = mincx(lmi_set,c,opts);


if (isempty(xopt))
    out.feas = 0;
    return
end
out.feas = 1;

out.W = dec2mat(lmi_set,xopt,W);
out.P = dec2mat(lmi_set,xopt,P);
%out.it = it;
end
