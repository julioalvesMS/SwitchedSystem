
%function affine_TimeVar_Ex1(kappa)
close all; clear; clc;

kappa = 10;

[F,g,E,G,H,T] = ExTheorical;
[A,b] = c2d_exact(F,g,T);    
x0= [1;1];

xe = [-9;0];

dsys.A = A;
dsys.b = b;
dsys.Q{1} = E{1}'*E{1};
dsys.x0 = [1;1];
Gamma = [1 0];

run system_specifications
circuit = buck(R, Ro, Co, L);
xe = [0; 30];
Gamma = [0 1];
T = Ts;
run load_circuit_sys

A = dsys.A;
b = dsys.b;
x0 = dsys.x0;
N = length(dsys.A);
n = length(dsys.A{1});

kappa = 13;
cand = find_cycles(dsys, xe, Gamma, kappa)

opt = 1;

cycle = find_limit_cycle(dsys, kappa, cand, opt);


xe_h = cycle.xe_h;
P = cycle.lyap.P;
ell = cycle.ell;

switch opt
    case 1
        sigma = @(i_,x_) sigma_h2(A, ell, xe_h, P, E, i_, x_(:,i_-1));
        wf = @(i_) 0;
    case 2
        sigma = @(i_,x_) sigma_h2(A, ell, xe_h, P, E, i_, x_(:,i_-1));
        wf = @(i_) 0;
        m=out.m;
        x0 = A{m}*xe_h{kappa}+b{m}+H{m};
    case 3
        rho = out.optval;
        sigma = @(i_,x_) sigma_hinf(A, ell, xe_h, P, E, H, G, rho, i_, x_(:,i_-1));
        wf = @(i_) 0.04*sin(2*pi*60*T*i_*(i_<.05/T));       
        x0 = xe_h{1};
end


K=floor(1/Ts); %[s]/T

x=zeros(n,K);
x(:,1) =x0;
xi=zeros(n,K);

for i=2:(K+1)
    s(i-1)= sigma(i,x);
    w(i-1) = wf(i);
    x(:,i) = A{s(i-1)}*x(:,i-1)+b{s(i-1)}+H{s(i-1)}*w(i-1) ;
    xi(:,i-1) = (x(:,i-1)-xe_h{mod(i-2,kappa)+1});
    ze(:,i-1) = E{s(i-1)}*xi(:,i-1)+G{s(i-1)}*w(i-1);
    v(i-1) = xi(:,i-1)'*P{mod(i-2,kappa)+1}*xi(:,i-1);
end

if(opt<2)   
        cost = trace(ze*ze');
elseif(opt<3)   
        cost = trace(ze*ze'+G{m}'*G{m});
else
        cost = trace(ze*ze')/trace(w*w')
end


f=figure(1);
f.Position(3:4) = [560 400];
plot(x(1,:),x(2,:),'LineWidth',.2,'Marker','.','MarkerSize',10)
hold on
for i=1:kappa
    scatter(xe_h{i}(1),xe_h{i}(2),'LineWidth',7,'Marker','*')
end
xlabel('$x_1[n]$','interpreter','latex','FontSize',12)
ylabel('$x_2[n]$','interpreter','latex','FontName','Times-Roman','FontSize',12)
%set(gca,'FontSize',15)

annotation(f,'arrow',[0.310714285714286 0.239285714285714],...
    [0.2065 0.2625],'LineStyle','none');


k= 0:(length(s)-1);
f=figure(2);
f.Position(3:4) = [560 240];
plot(k, xi','LineWidth',.2,'Marker','.','MarkerSize',10)
xlabel('$n$','interpreter','latex','FontName','Times-Roman','FontSize',12)
ylabel('$\xi[n]$','interpreter','latex','FontName','Times-Roman','FontSize',12)
ylim([-1 15])
xlim([0 70])
%set(gca,'FontSize',15)

f=figure(3);
f.Position(3:4) = [560 240];
stairs(k,s,'LineWidth',.2,'Marker','.','MarkerSize',10)
yticks([1 2])
ylim([.8 2.2])
xlim([0 70])
xlabel('$n$','interpreter','latex','FontName','Times-Roman','FontSize',12)
ylabel('$\sigma[n]$','interpreter','latex','FontName','Times-Roman','FontSize',12)
%axis([0 K 0.95 2.05])
%set(gca,'FontSize',15)

%end


%%
function [s] = sigma_h2(A,ell,xe,P,E,n,x) 
    kappa= length(P);
    k =mod(n-2,kappa)+1;
    km=mod(n-1,kappa)+1;
    for i=1:length(A)
        val(i) = [x-xe{k};1]'*([A{i}'*P{km}*A{i}-P{k}+E{i}'*E{i}  A{i}'*P{km}*ell{i,k};  ell{i,k}'*P{km}*A{i} ell{i,k}'*P{km}*ell{i,k}]*[x-xe{k};1]);
    end
    [~,s] = min(val);
end
%%
function [s] = sigma_hinf(A,ell,xe,P,E,H,G,rho,n,x) 
    kappa= length(P);
    k =mod(n-2,kappa)+1;
    km=mod(n-1,kappa)+1;
    for i=1:length(A)
        Lcal = [A{i}'*P{km}*A{i}-P{k}      A{i}'*P{km}*ell{i,k}; nm
                ell{i,k}'*P{km}*A{i}        ell{i,k}'*P{km}*ell{i,k}];
        Dcal = [A{i}'*P{km}*H{i}+E{i}'*G{i};
                ell{i,k}'*P{km}*H{i}];
        Sigmacal = H{i}'*P{km}*H{i}+G{i}'*G{i}-rho*eye(size(H{1},2));
        Mcal=Dcal/Sigmacal*Dcal';
        val(i) = [x-xe{k};1]'*(Lcal-Mcal)*[x-xe{k};1]-...
        (x-xe{k})'*E{i}'*E{i}*(x-xe{k});
    end
    [~,s] = min(val);
end

%% 
function [F,g,E,G,H,T] = ExTheorical


F{1} = [-4 3;-3 2.5];
F{2} = [4 -1;1 -2];

g{1} = [0;-2];
g{2} = [0;8];

T = 1e-1;

E{1} = eye(2);
E{2} = E{1};

G{1} = [0;0];
G{2} = G{1};


H{1} = [0;0];
H{2} = [0;0];
end
