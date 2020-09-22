%%Ex1_main 
% This program runs Example 1 optmization problems and output table 1 data
% mod: theta based
function Ex1_main
%System data
ExIEEEdata
%ExECCdata
%Discretization
[A,b] = c2d_exact(F,g,T);
Al = convc(A,lambda);
bl = convc(b,lambda);
%System state-space deplacement
xe = (eye(n) - Al)\bl;
for i =1:N
    ell{i} = (A{i}-eye(n))*xe+b{i};
end
% 
Theta = [0 1e5; 1 0];
Theta = Theta/max(eig(Theta))
% %Solving Theorem 1 (SF)
 [out_T1]= opt_theorem1(A,ell,Theta);

 
%  
% % %Solving Corollary 1 (SF)
%  [out_C1]= opt_theorem1(A,ell,MpiC1,ones(1,N));
% % 
% % %Solving Theorem 2 (SF) for both T1 and C1
% % [out_T2_T1] = opt_theorem2(A,ell,MpiT1,out_T1.P);
% % [out_cvx_inv_T1] = min_union(out_T1.P,out_T2_T1.r);
% % 
% % [out_T2_C1] = opt_theorem2(A,ell,MpiC1,out_C1.P);
% % [out_cvx_inv_C1] = min_union(out_C1.P,out_T2_C1.r);
% 
% [out_T2_T1] = opt_theorem2_v2(A,ell,MpiT1,out_T1.P);
% [out_cvx_inv_T1] = min_union(out_T1.P,out_T2_T1.r);
% 
% [out_T2_C1] = opt_theorem2_v2(A,ell,MpiC1,out_C1.P);
% [out_cvx_inv_C1] = min_union(out_C1.P,out_T2_C1.r);
% 
% % 
% % %Solving IEEE T1,T2 (SF)
% [out_IEEE_T1] = opt_IEEEtheorem1(A,ell,lambda);
% [out_IEEE_T2] = opt_IEEEtheorem2(out_IEEE_T1.h,out_IEEE_T1.P,out_IEEE_T1.W);
% 
% %Solving IEEE T3 (SF)
% [out_IEEE_T3] = opt_IEEEtheorem3(A,ell,lambda);
% 
% 
% 
% %Solving ECC  (SF)
% [out_ECC]= opt_ECC(A,ell,lambda);
% 
% %Solving Theorem 3 (OF)
% [out_OF]= opt_OF(A,ell,C,MpiOF,betasOF);
% 
% 
% %Solving Corollary 2 (OF)
% [out_OFC]= opt_OFC(A,ell,C,MpiOF) ;

%Outputs Table 1
keyboard
     
Table = [nvol(out_T1.W,n)                                       nvol(out_cvx_inv_T1.R,n);
         nvol(out_C1.W,n)                                       nvol(out_cvx_inv_C1.R,n);
         nvol(out_IEEE_T1.W,n),                                 nvol(out_IEEE_T1.P/out_IEEE_T2.r,n);
         nvol(out_IEEE_T3.W/convc(out_IEEE_T3.rho,lambda),n)    nvol(out_IEEE_T3.P,n)
         nvol(out_ECC.W,n)                                      NaN;
         nvol(inv(out_OF.W),n)                                  NaN 
         nvol(inv(out_OFC.W)*out_OFC.beta^2,n)                                  NaN]%/nvol(out_T1.W,n)*100
%



figure(1)
xc_i  = -out_union.M\out_union.a;
[x,y,z]=plot_ellipsoid(out_union.M,xc_i,1);


keyboard
end


%% opt_theorem1
function [out]= opt_theorem1(A,ell,Theta) 
N = length(A);
n = length(A{1});
setlmis([])

%Variables
W = lmivar(1,[n 1]);
h = lmivar(2,[n 1]);
for i =1:N
    S{i} = lmivar(1,[n 1]);
end
    
for i = 1:N
    %Sprocedure
    clmi = newlmi;
    %   Col1
    lmiterm([-clmi,1,1,S{i}],1,1);
    for j=1:N
        lmiterm([-clmi,3,1,S{j}],Theta(j,i),A{i});
    end
    lmiterm([-clmi,4,1,W],-1,1);
    %   Col2
    lmiterm([-clmi,2,2,0],1);
    for j=1:N
        lmiterm([-clmi,3,2,S{j}],Theta(j,i),ell{i});
    end
    lmiterm([-clmi,4,2,h],-1,1);
    %   Col3
    for j=1:N
        lmiterm([-clmi,3,3,S{j}],Theta(j,i),1);
    end
    %   Col4
    lmiterm([-clmi,4,4,W],1,1);
end

%Solution
lmi_set = getlmis;
opts = [1e-4,2000,-1,10,1];
[optval,xopt,it] = frank_wolfe(lmi_set,opts,W);

if (isempty(xopt)) % No Feas ?
    out.feas = 0;
    return
end
[V,~] = eig(Theta');
out.feas = 1; % Feas
out.W = dec2mat(lmi_set,xopt,W);
out.h = dec2mat(lmi_set,xopt,h);
for i=1:N
    out.P{i} = dec2mat(lmi_set,xopt,S{i})/abs(V(i,1));
end
out.it = it;
end

%% opt_theorem1_v2
% takes into account sigma=i in volume minimization
function [out]= opt_theorem1_v2(A,ell,Mpi,betas) 
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
[optval,xopt,it] = frank_wolfe(lmi_set,opts,W);

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
out.it = it;
end
%% opt_theorem2
function [out] = opt_theorem2(A,ell,Mpi,P,idx)
setlmis([])
N = length(P);
span = 1:N;
if(exist('idx','var'))
    span =idx;
end
%Variables
r = lmivar(1,[1 0]);
for i = span
    alphas(i)  = lmivar(1,[1 0]);
end
 
for i = span
    %Sprocedure
    Ppi = convc(P,Mpi(:,i));
    clmi = newlmi;
    %   Col1
    lmiterm([-clmi,1,1,alphas(i)],P{i},1);
    lmiterm([-clmi,1,1,0],-P{i});
    lmiterm([-clmi,3,1,alphas(i)],1,Ppi*A{i});
    
    %   Col2
    lmiterm([-clmi,2,2,r],1,1);
    lmiterm([-clmi,3,2,alphas(i)],1,Ppi*ell{i});
    %   Col3
    lmiterm([-clmi,3,3,alphas(i)],1,Ppi);
end

lmisys = getlmis;
%===============

options = [1e-4,2000,-1,200,1];
np = decnbr(lmisys);

c = zeros(np,1);
c(decinfo(lmisys,r)) = 1;

[copt,xopt] = mincx(lmisys,c,options);

if (isempty(copt))
    out.feas = 0;
else
    out.feas =1;
    out.r =  dec2mat(lmisys,xopt,r);
    for i = span
        out.alpha(i) = dec2mat(lmisys,xopt,alphas(i));
    end
end
end

%% min_union
function [out] = min_union(P,r)
setlmis([])
N = length(P);
if(length(r)==1)
   r= r*ones(1,N); 
end
n  =length(P{1});
R = lmivar(1,[n 1]);
d = lmivar(2,[n 1]);
for i=1:N
    theta(i)  = lmivar(1,[1 0]);
end


for i =1:N
    clmi = newlmi;
    %   Col1
    lmiterm([-clmi,1,1,R],-1,1);
    lmiterm([-clmi,1,1,theta(i)],P{i},1);
    
    lmiterm([-clmi,2,1,-d],1,1);
    %   Col2
    lmiterm([-clmi,2,2,0],1);
    lmiterm([-clmi,2,2,theta(i)],-r(i),1);
    lmiterm([-clmi,3,2,d],1,1);
    %   Col3
    lmiterm([-clmi,3,3,R],1,1);

    
    clmi = newlmi;
    lmiterm([-clmi,1,1,theta(i)],1,1);
end


lmisys = getlmis;
options = [1e-4,2000,-1,200,1];

[copt,xopt,it] = frank_wolfe(lmisys,options,R);

if (isempty(copt))
    out.feas = 0;
else
    out.feas =1;
    out.R =  dec2mat(lmisys,xopt,R);
    out.d =  dec2mat(lmisys,xopt,d);
    for i =1:N
        out.theta(i) = dec2mat(lmisys,xopt,theta(i));
    end
    out.it = it;
end
end


%% frank_wolfe
function [optval,xopt,it] = frank_wolfe(lmi_set,opts,Rvar)
epsilon = 1e-3;

it = 1;
n = length(decinfo(lmi_set, Rvar));
Rx = eye(n);


[tmin,xopt] = feasp(lmi_set,opts);

if tmin>0
    xopt = [];
    optval = [];
    return;
else
    
    Rx  = dec2mat(lmi_set,xopt,Rvar);
    outx = xopt;
    
    while 1
        n_dec = decnbr(lmi_set);
        c = zeros(n_dec,1);
        for i=1:n_dec
            Rn =  defcx(lmi_set,i,Rvar);  
            c(i) = -trace(Rx\Rn );
        end
        

        [~,xopt] = mincx(lmi_set,c,opts);
        out = xopt;
        if( isempty(xopt))
            keyboard
        end
        R = dec2mat(lmi_set,xopt,Rvar);
        cond = -trace((Rx)\(R-Rx)) ;
     %   fprintf('%4d: %f\n',it, cond)  %show converg.
        if cond < -epsilon
            vet = @(x_) -log(det(x_*R+(1-x_)*Rx));
            
            [alphax,optval]= unimin(vet);
            Rx =   alphax*R   + (1-alphax)*Rx;
     
            outx = alphax*out + (1-alphax)*outx;
            
        else
            
            xopt =outx;
            break
        end
        it = it+1;
    end
end
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


[copt,xopt,it] = frank_wolfe(lmisys,options,W);


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

%% opt_IEEEtheorem2
function [out] = opt_IEEEtheorem2(h,P,W)
setlmis([])

%Variables
r = lmivar(1,[1 0]);
betav = lmivar(1,[1 0]);
 
%Sprocedure
clmi = newlmi;
%   Lin1
lmiterm([-clmi,1,1,r],1,1);
lmiterm([-clmi,1,1,betav],-1,1);
lmiterm([-clmi,1,2,0],h');

%   Lin2
lmiterm([-clmi,2,2,0],P);
lmiterm([-clmi,2,3,0],P);
%   Lin3
lmiterm([-clmi,3,3,betav],1,W);


lmisys = getlmis;
%===============

options = [1e-4,2000,-1,200,1];
np = decnbr(lmisys);

c = zeros(np,1);
c(decinfo(lmisys,r)) = 1;

[copt,xopt] = mincx(lmisys,c,options);

if (isempty(copt))
    out.feas = 0;
else
    out.feas =1;
    out.r =  dec2mat(lmisys,xopt,r);
    out.beta = dec2mat(lmisys,xopt,betav);
end
end


%%  opt_IEEEtheorem3
function [out] = opt_IEEEtheorem3(A,ell,lambin)

betaInt = [1e-4 1e10];

[x,fval]= unimin(@(beta) AUX_opt_IEEEtheorem3(A,ell,lambin,beta,0),1e-2,0,betaInt);

 [out] = AUX_opt_IEEEtheorem3(A,ell,lambin,x,1);
end
%%  AUX_opt_IEEEtheorem3
function [out] = AUX_opt_IEEEtheorem3(A,ell,lambin,beta,outputType)
n = size(A{1},1);
N = size(A,2);
setlmis([])

%Variables
P = lmivar(1,[n 1]);
W = lmivar(1,[n 1]);
               

Al = convc(A,lambin);
clmi = newlmi; 
for i = 1:N
    lmiterm([-clmi 1 1 P],lambin(i)*beta*ell{i}',-ell{i});                  % LMI #1: -B*ell{i}'*P*ell{i} 
    lmiterm([-clmi 2 1 P],lambin(i)*((A{i}/(eye(n)-Al)))',(-ell{i}));       % LMI #1: (-ell{i}'*P*(A{i}/(eye(n)-Al)))'
end
lmiterm([-clmi 1 1 0],1);                                                   % LMI #1: 1
lmiterm([-clmi 2 2 P],1,1);                                                   % LMI #1: P
lmiterm([-clmi 3 2 P],1,1);                                                % LMI #1: P
lmiterm([-clmi 3 3 W],beta,1);                                              % LMI #1: B*W 


clmi = newlmi; 
for i=1:N
    lmiterm([clmi,1,1,P],lambin(i)*A{i}',A{i});
end
lmiterm([clmi,1,1,P],-1,1);
lmiterm([clmi,1,1,W],1,1);


lmisys = getlmis;
%===============



options = [1e-4,2000,-1,200,1];


[copt,xopt,it] = frank_wolfe(lmisys,options,P);


if (isempty(copt))
    if(outputType==1)
        out.feas=0;
    else
        out = inf;
    end
else
    if(outputType==1)
        out.feas=1;
        out.it= it;
        out.W = dec2mat(lmisys,xopt,W);
        out.P = dec2mat(lmisys,xopt,P);
        out.beta = beta;
        for i =1:N
            out.c{i} = A{i}'*out.P*ell{i};
            out.rho{i} = ell{i}'*out.P*ell{i};
        end
        out.h = (eye(n)-convc(A,lambin)')\convc(out.c,lambin);
    else
        out = -log(det(dec2mat(lmisys,xopt,P)));
    end
end
end


%% opt_theorem2
% takes into account \xi \in \min(v(xi))
function [out] = opt_theorem2_v2(A,ell,Mpi,P,idx)
setlmis([])
N = length(P);
span = 1:N;
if(exist('idx','var'))
    span =idx;
end
%Variables
r = lmivar(1,[1 0]);
for i = span
    alphas(i)  = lmivar(1,[1 0]);
    for j=span
        if(j==i)
            continue
        end
        phis(i,j)  = lmivar(1,[1 0]);
    end
end
 
for i = span
    %Sprocedure
    Ppi = convc(P,Mpi(:,i));
    clmi = newlmi;
    %   Col1
    lmiterm([-clmi,1,1,alphas(i)],P{i},1);
    for j=span
        if(i~=j)
            lmiterm([-clmi,1,1,phis(i,j)],P{i},1); 
            lmiterm([-clmi,1,1,phis(i,j)],-P{j},1); 
        end
    end
    lmiterm([-clmi,1,1,0],-P{i});
    lmiterm([-clmi,3,1,alphas(i)],1,Ppi*A{i});
    
    %   Col2
    lmiterm([-clmi,2,2,r],1,1);
    lmiterm([-clmi,3,2,alphas(i)],1,Ppi*ell{i});
    %   Col3
    lmiterm([-clmi,3,3,alphas(i)],1,Ppi);
    
    
    for j=span
        if(i==j)
            continue;
        end
        clmi = newlmi;

        lmiterm([-clmi,1,1,phis(i,j)],1,1);

    end
end
lmisys = getlmis;
%===============

options = [1e-4,2000,-1,200,1];
np = decnbr(lmisys);

c = zeros(np,1);
c(decinfo(lmisys,r)) = 1;

[copt,xopt] = mincx(lmisys,c,options);

if (isempty(copt))
    out.feas = 0;
else
    out.feas =1;
    out.r =  dec2mat(lmisys,xopt,r);
    for i = span
        out.alpha(i) = dec2mat(lmisys,xopt,alphas(i));
        for j=span
            if (i==j)
                continue
            end
            out.phis(i,j) = dec2mat(lmisys,xopt,phis(i,j));
        end
    end
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



%Solucao
lmi_set = getlmis;

opts = [1e-4,2000,-1,10,1];


[optval,xopt,it] = frank_wolfe(lmi_set,opts,W);


if (isempty(xopt))
    out.feas = 0;
    return
end
out.feas = 1;

out.W = dec2mat(lmi_set,xopt,W);
out.P = dec2mat(lmi_set,xopt,P);
out.it = it;
end

%% opt_OF
function [out]= opt_OF(A,ell,C,Mpi,beta) 
N = length(A);
n = length(A{1});
m = size(C{1},1);
setlmis([])

%Variaveis
Y = lmivar(1,[n 1]);
W = lmivar(1,[n 1]);

g = lmivar(2,[n 1]);
for i =1:N
    nu{i} = lmivar(2,[n 1]);
    L{i} = lmivar(2,[n m]);
    
    X{i} = lmivar(1,[n 1]);
    M{i} = lmivar(2,[n n]);
    J{i} = lmivar(2,[n n]);
    for j=1:N
        T{i,j} = lmivar(1,[n 1]);
    end
end


for i = 1:N
    %LM
    clmi = newlmi;
    lmiterm([-clmi,1,1,Y],beta(i),1);
    lmiterm([-clmi,2,1,0],beta(i));
    lmiterm([-clmi,4,1,Y],beta(i),A{i});
    lmiterm([-clmi,4,1,L{i}],beta(i),C{i});
    lmiterm([-clmi,5,1,0],beta(i)*A{i});
    lmiterm([-clmi,6,1,0],1);
    
    
    lmiterm([-clmi,2,2,X{i}],beta(i),1);
    lmiterm([-clmi,4,2,M{i}],beta(i),1);
    lmiterm([-clmi,5,2,X{i}],A{i},beta(i));
    lmiterm([-clmi,6,2,X{i}],1,1);
    
    
    lmiterm([-clmi,3,3,0],1);
    lmiterm([-clmi,4,3,Y],beta(i),ell{i});
    lmiterm([-clmi,4,3,nu{i}],beta(i),1);
    lmiterm([-clmi,5,3,0],beta(i)*ell{i});
    lmiterm([-clmi,6,3,g],-1,1);
    
    lmiterm([-clmi,4,4,Y],beta(i),1);
    lmiterm([-clmi,5,4,0],beta(i));

      
    lmiterm([-clmi,5,5,J{i}],beta(i),1,'s');
    for j=1:N
        lmiterm([-clmi,5,5,T{i,j}],-Mpi(j,i),beta(i));
    end
    
    lmiterm([-clmi,6,6,W],1,1);
  
    
    for j=1:N
        clmi = newlmi;
        lmiterm([-clmi,1,1,T{i,j}],1,1);
        lmiterm([-clmi,2,1,J{i}],1,1);
        lmiterm([-clmi,3,1,0],1);
        
        lmiterm([-clmi,2,2,X{j}],1,1);
        lmiterm([-clmi,3,2,0],1);
        
        lmiterm([-clmi,3,3,Y],1,1);
    
    end
    
      
    %No degenerated Yci
%     clmi = newlmi;
%     lmiterm([-clmi,1,1,Y],1,1);
%     %lmiterm([-clmi,1,1,0],-1e-2);
%     lmiterm([-clmi,2,1,0],1);
%     lmiterm([-clmi,2,2,X{i}],1,1);

end

     

%Solution
lmi_set = getlmis;

opts = [1e-4,2000,-1,10,1];

n_dec = decnbr(lmi_set);
c = zeros(n_dec,1);
for i=1:n_dec
    Wn =  defcx(lmi_set,i,W);
    c(i) = trace(Wn);
end

it =1;
[~,xopt] = mincx(lmi_set,c,opts);
%[optval,xopt,it] = frank_wolfe(lmi_set,opts,[R1, R2,R3]);


if (isempty(xopt))
    out.feas = 0;
    return
end
out.feas = 1;






out.Y = dec2mat(lmi_set,xopt,Y);
out.W = dec2mat(lmi_set,xopt,W);
out.g = dec2mat(lmi_set,xopt,g);
for i=1:N
    out.nu{i} = dec2mat(lmi_set,xopt,nu{i});
    out.L{i} = dec2mat(lmi_set,xopt,L{i});
    out.X{i} = dec2mat(lmi_set,xopt,X{i});
    out.M{i} = dec2mat(lmi_set,xopt,M{i});
    out.J{i} = dec2mat(lmi_set,xopt,J{i});
    for j=1:N
        out.T{i,j} = dec2mat(lmi_set,xopt,T{i,j});
    end
end

out.it = it;
end

%% opt_OFC
function [out]= opt_OFC(A,ell,C,Mpi) 
N = length(A);
n = length(A{1});
m = size(C{1},1);
setlmis([])

%Variaveis
beta = lmivar(1,[1 1]);
Y = lmivar(1,[n 1]);
W = lmivar(1,[n 1]);

g = lmivar(2,[n 1]);
for i =1:N
    nu{i} = lmivar(2,[n 1]);
    L{i} = lmivar(2,[n m]);
    
    X{i} = lmivar(1,[n 1]);
    M{i} = lmivar(2,[n n]);
    J{i} = lmivar(2,[n n]);
    for j=1:N
        T{i,j} = lmivar(1,[n 1]);
    end
end


for i = 1:N
    %LM
    clmi = newlmi;
    lmiterm([-clmi,1,1,Y],1,1);
    lmiterm([-clmi,2,1,beta],1,1);
    lmiterm([-clmi,4,1,Y],1,A{i});
    lmiterm([-clmi,4,1,L{i}],1,C{i});
    lmiterm([-clmi,5,1,beta],1,A{i});
    lmiterm([-clmi,6,1,beta],1,1);
    
    
    lmiterm([-clmi,2,2,X{i}],1,1);
    lmiterm([-clmi,4,2,M{i}],1,1);
    lmiterm([-clmi,5,2,X{i}],A{i},1);
    lmiterm([-clmi,6,2,X{i}],1,1);
    
    
    lmiterm([-clmi,3,3,0],1);
    lmiterm([-clmi,4,3,Y],1,ell{i});
    lmiterm([-clmi,4,3,nu{i}],1,1);
    lmiterm([-clmi,5,3,beta],1,ell{i});
    lmiterm([-clmi,6,3,g],-1,1);
    
    lmiterm([-clmi,4,4,Y],1,1);
    lmiterm([-clmi,5,4,beta],1,1);

      
    lmiterm([-clmi,5,5,J{i}],1,1,'s');
    for j=1:N
        lmiterm([-clmi,5,5,T{i,j}],-Mpi(j,i),1);
    end
    
    lmiterm([-clmi,6,6,W],1,1);
  
    
    for j=1:N
        clmi = newlmi;
        lmiterm([-clmi,1,1,T{i,j}],1,1);
        lmiterm([-clmi,2,1,J{i}],1,1);
        lmiterm([-clmi,3,1,beta],1,1);
        
        lmiterm([-clmi,2,2,X{j}],1,1);
        lmiterm([-clmi,3,2,beta],1,1);
        
        lmiterm([-clmi,3,3,Y],1,1);
    
    end
    
      
    %No degenerated Yci
%     clmi = newlmi;
%     lmiterm([-clmi,1,1,Y],1,1);
%     %lmiterm([-clmi,1,1,0],-1e-2);
%     lmiterm([-clmi,2,1,0],1);
%     lmiterm([-clmi,2,2,X{i}],1,1);

end

    
     clmi = newlmi; %beta>1e-3
     lmiterm([-clmi,1,1,beta],1,1);
     lmiterm([-clmi,1,1,0],-1e-3);
     

%Solution
lmi_set = getlmis;

opts = [1e-4,2000,-1,10,1];

n_dec = decnbr(lmi_set);
c = zeros(n_dec,1);
for i=1:n_dec
    Wn =  defcx(lmi_set,i,W);
    c(i) = trace(Wn);
end

it =1;
[~,xopt] = mincx(lmi_set,c,opts);
%[optval,xopt,it] = frank_wolfe(lmi_set,opts,[R1, R2,R3]);


if (isempty(xopt))
    out.feas = 0;
    return
end
out.feas = 1;





out.beta = dec2mat(lmi_set,xopt,beta);

out.Y = dec2mat(lmi_set,xopt,Y);
out.W = dec2mat(lmi_set,xopt,W);
out.g = dec2mat(lmi_set,xopt,g);
for i=1:N
    out.nu{i} = dec2mat(lmi_set,xopt,nu{i});
    out.L{i} = dec2mat(lmi_set,xopt,L{i});
    out.X{i} = dec2mat(lmi_set,xopt,X{i});
    out.M{i} = dec2mat(lmi_set,xopt,M{i});
    out.J{i} = dec2mat(lmi_set,xopt,J{i});
    for j=1:N
        out.T{i,j} = dec2mat(lmi_set,xopt,T{i,j});
    end
end

out.it = it;
end




