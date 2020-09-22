function Ex1_simulation_main
    close all
    %System data
    %Ex2dNovo
    ExLMdata
    %Discretization
    [A,b] = c2d_exact(F,g,T);
    Al = convc(A,lambda);
    bl = convc(b,lambda);
    %System state-space deplacement
    
    xe = (eye(n) - Al)\bl
    for i =1:N
        ell{i} = (A{i}-eye(n))*xe+b{i};
    end
    
    %Solving Theorem 1 and 2
    [out_T1]= opt_theorem1(A,ell,MpiT1,betas);
    [out_C1]= opt_theorem1(A,ell,MpiC1,[1 1]);
    [out_T2_T1] = opt_theorem2_v2(A,ell,MpiT1,out_T1.P);
    
    keyboard
    %%
     close all
    %---Graphical Output - Figure 1
    P = out_T1.P;
    W = out_T1.W;
    h = out_T1.h;

    for i=1:N
        Ppi{i} = P{1}*0;
        for j=1:N
            Ppi{i} = Ppi{i} + P{j}*MpiT1(j,i);
        end
        Q{i}    = -A{i}'*Ppi{i}*A{i}+P{i};
        c{i}    = A{i}'*Ppi{i}*ell{i};
        rho{i}  = ell{i}'*Ppi{i}*ell{i};

    end

    
    
    Tfin = 50;
    Ntime = floor(Tfin/T);
    tspan = linspace(0,Tfin-T,Ntime);
    xi0 = 20*[+1 -1 -1]';
    [xi,sigma] = time_simulation(A, ell,P,Ntime,xi0);
    for i=1:length(xi)
        lyap(i) =  (xi(:,i))'*P{sigma(i)} *(xi(:,i));
    end
    %% option 1
    close all
    f=figure(1);
    subplot(2,1,1);
    h= stairs(tspan,lyap);
    set(h.get('Parent'),'Yscale', 'log');
    ylabel('$v(\xi)$');
    hold on
    h = plot(tspan,tspan*0+out_T2_T1.r,'k--')
 
    leg = legend(h,['$r_\ast$'])
    leg.Position=[0.633304380415249 0.85 0.255057905836168 0.0770000000000001];
   % leg= legend('$\xi_1$', '$\xi_2$', '$\xi_3$');
    leg.Interpreter = 'latex';
    legend('boxoff')
    set( subplot(2,1,1),'YMinorTick','on','YScale','log','YTick',[1 100000],...
    'YTickLabel',{'10^{0}','10^5'});

    xlim([0 20])
    subplot(2,1,2);
    stairs(tspan,sigma);
    xlim([0 20])
    ylim([0.5 2.5]);
    ylabel('$\sigma$');
    xlabel('Time [s]');
    yticks([1 2])
    set(f,'Position', [200 200 1600/2 1100]*0.3);
    
    
    f = figure(2);
    set(f,'Position', [200 200 1600/2 1100]*0.3);
     h_3d=plot3(xi(1,:),xi(2,:),xi(3,:),'Marker','o','LineWidth',2)
     axes1 = get(h_3d, 'Parent');
    hold on
    for i =1:N
        [x,y,z]=plot_ellipsoid(out_T1.P{i},xe*0,out_T2_T1.r);
        mesh(x,y,z,'FaceLighting','none',...
    'EdgeLighting','flat',...
    'FaceAlpha',0.3,...
    'FaceColor',[1 0 0],...
    'EdgeColor','none');
    end
    
xlabel('$\xi_1$');

zlabel('$\xi_3$');

ylabel('$\xi_2$');

    xlim([-20 40])
    view(axes1,[47.7 51.6000000000001]);
box(axes1,'on');
grid(axes1,'on');


%% option 2
close all
    f=figure(3);
    ax1=subplot(2,1,1);
    h= plot(tspan,lyap,'Marker','.','MarkerSize',14);
    set(h.get('Parent'),'Yscale', 'log');
    ylabel('vxin');
    hold on
    h = plot(tspan,tspan*0+out_T2_T1.r,'k--')
 
    leg = legend(h,['ra'])
    leg.Position=[0.6922 0.8381 0.2551 0.0770];
   % leg= legend('$\xi_1$', '$\xi_2$', '$\xi_3$');
    %leg.Interpreter = 'latex';
    legend('boxoff')
    set( subplot(2,1,1),'YMinorTick','on','YScale','log','YTick',[1 100000],...
    'YTickLabel',{'10^{0}','10^5'});

    xlim([0 20])
    ax2=subplot(2,1,2);
    stairs(tspan,sigma,'Marker','.','MarkerSize',14);
    xlim([0 20])
    ylim([0.5 2.5]);
    ylabel('sig');
    xlabel('n');
    yticks([1 2])
   % set(f,'Position', [200 200 1600/2 1100]*0.3);
    f3 = 'ex2_lyap1.pdf';
opts.scriptsize=1;
opts.frags= { {'vxin','$v(\xi[n])$'};
              {'sig','$\sigma[n]$'};
              {'n','$n$'}}
saveFigLatex(f3,opts)
    
    f = figure(4);
 %   set(f,'Position', [200 200 1600/2 1100]*0.3);
     h_3d=plot3(xi(1,:),xi(2,:),xi(3,:),'Marker','.','MarkerSize',14,'LineWidth',2)
     axes1 = get(h_3d, 'Parent');
     ax3= axes1
    hold on
    for i =1:N
        [x,y,z]=plot_ellipsoid(out_T1.P{i},xe*0,out_T2_T1.r);
        mesh(x,y,z,'FaceLighting','none',...
    'EdgeLighting','flat',...
    'FaceAlpha',0.3,...
    'FaceColor',[1 0 0],...
    'EdgeColor','none');
    end
    xlim([-20 40])
    view(axes1,[47.7 51.6000000000001]);
box(axes1,'on');
grid(axes1,'on');
ax1.FontSize = 12;
ax2.FontSize = 12;
ax3.FontSize = 12;
       
    
xlabel('$\xi_1$','Interpreter','latex','FontSize',16);

ylabel('$\xi_2$','Interpreter','latex','FontSize',16);

zlabel('$\xi_3$','Interpreter','latex','FontSize',16);

%%
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
    lmiterm([-clmi,4,1,W],1,1);
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
%% time_simulation
function [xi,sigma] = time_simulation(A, ell ,P ,K,xi0) 
    N = length(A);
    n = length(xi0);
    xi = zeros(n,K);
    sigma = zeros(1,K)+1;
    xi(:,1)=xi0;
    
    for k = 2:K
      %  fprintf('%.2f  \n',k/K*100) %show progress
        for i=1:N
            val(i) = xi(:,k-1)'*P{i}*xi(:,k-1);
       end            
        [~,i] = min(val);
        sigma(k-1)=i;
        xi(:,k) = A{i}*xi(:,k-1)+ell{i};    
    end
    for i=1:N
            val(i) = xi(:,K)'*P{i}*xi(:,K);
    end            
    [~,i] = min(val);
    sigma(K)=i;




end


%% plot_region
function [x,y] = plot_region(W,hw,rw, Delt)
theta_grid = linspace(0,2*pi,2000);
% el1
P =  W{1}/rw{1};
r = (eig(P)).^(-1/2);
[rot, ~] = eig(P);


phi = atan2(rot(2,1), rot(1,1));
ellipse_x_r  = r(1)*cos( theta_grid );
ellipse_y_r  = r(2)*sin( theta_grid );
R = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];
r_ellipse = [ellipse_x_r;ellipse_y_r]' * R;
x1=r_ellipse(:,1)+hw{1}(1);
y1=r_ellipse(:,2)+hw{1}(2);

% el2
P =  W{2}/rw{2};
r = (eig(P)).^(-1/2);
[rot, ~] = eig(P);


phi = atan2(rot(2,1), rot(1,1));
ellipse_x_r  = r(1)*cos( theta_grid );
ellipse_y_r  = r(2)*sin( theta_grid );
R = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];
r_ellipse = [ellipse_x_r;ellipse_y_r]' * R;
x2=r_ellipse(:,1)+hw{2}(1);
y2=r_ellipse(:,2)+hw{2}(2);

for i = 1:length(x1)
    if([x1(i) y1(i)]*Delt*[x1(i);y1(i)]>0)
        x1(i)=NaN;
        y1(i)=NaN;
    end
    if([x2(i) y2(i)]*Delt*[x2(i);y2(i)]<0)
        x2(i)=NaN;
        y2(i)=NaN;
    end
end
%x = [x1(~isnan(x1)); x2(~isnan(x2))];
%y = [y1(~isnan(y1)); y2(~isnan(y2))];

x1 = splitNaN(x1');
y1 = splitNaN(y1');
x2 = splitNaN(x2');
y2 = splitNaN(y2');


for i = 1:(length(x1)-1)
    val_aux = 0*[1:length(x2)];
    for j = 1:length(x2)
        aux = (cross([x1{i}(end) y1{i}(end) 0],[x2{j}(1) y2{j}(1) 0]));
        if(aux(3)<0)
            aux(3) = inf;
        end
        val_aux(j)  = aux(3);
    end
    [~, idx]=min(val_aux);
    x1{i}(end+1) = x2{idx}(1);
    y1{i}(end+1) = y2{idx}(1);
end


for i = 1:(length(x2)-1)
    val_aux = 0*[1:length(x1)];
    for j = 1:length(x1)
        aux = (cross([x2{i}(end) y2{i}(end) 0],[x1{j}(1) y1{j}(1) 0]));
        if(aux(3)<0)
            aux(3) = inf;
        end
        val_aux(j)  = aux(3);
    end
    [~, idx]=min(val_aux);
    x2{i}(end+1) = x1{idx}(1);
    y2{i}(end+1) = y1{idx}(1);
end


x = [];
y = [];
for i = 1:length(x1)
    x = [x NaN x1{i}];
    y = [y NaN y1{i}];
end
for i = 1:length(x2)
    x = [x NaN x2{i}];
    y = [y NaN y2{i}];
end
end

%% plot_region2
function [x,y] = plot_region2(W,hw,rw)
theta_grid = linspace(0,2*pi,2000);
% el1
P =  W{1}/rw{1};
r = (eig(P)).^(-1/2);
[rot, ~] = eig(P);

phi = atan2(rot(2,1), rot(1,1));
ellipse_x_r  = r(1)*cos( theta_grid );
ellipse_y_r  = r(2)*sin( theta_grid );
R = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];
r_ellipse = [ellipse_x_r;ellipse_y_r]' * R;
x1=r_ellipse(:,1)+hw{1}(1);
y1=r_ellipse(:,2)+hw{1}(2);

% el2
P =  W{2}/rw{2};
r = (eig(P)).^(-1/2);
[rot, ~] = eig(P);


phi = atan2(rot(2,1), rot(1,1));
ellipse_x_r  = r(1)*cos( theta_grid );
ellipse_y_r  = r(2)*sin( theta_grid );
R = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];
r_ellipse = [ellipse_x_r;ellipse_y_r]' * R;
x2=r_ellipse(:,1)+hw{2}(1);
y2=r_ellipse(:,2)+hw{2}(2);

funcpert = @(xx,n) (xx-hw{n})'*W{n}*(xx-hw{n})<rw{n};

for i = 1:length(x1)
    if(funcpert([x1(i) y1(i)]',2))
        x1(i)=NaN;
        y1(i)=NaN;
    end
    if(funcpert([x2(i) y2(i)]',1))
        x2(i)=NaN;
        y2(i)=NaN;
    end
end
%x = [x1(~isnan(x1)); x2(~isnan(x2))];
%y = [y1(~isnan(y1)); y2(~isnan(y2))];
% 

x = [x1' NaN x2'];
y = [y1' NaN y2'];
% x1 = splitNaN(x1');
% y1 = splitNaN(y1');
% x2 = splitNaN(x2');
% y2 = splitNaN(y2');
% 
% 
% for i = 1:(length(x1)-1)
%     val_aux = 0*[1:length(x2)];
%     for j = 1:length(x2)
%         aux = (cross([x1{i}(end) y1{i}(end) 0],[x2{j}(1) y2{j}(1) 0]));
%         if(aux(3)<0)
%             aux(3) = inf;
%         end
%         val_aux(j)  = aux(3);
%     end
%     [~, idx]=min(val_aux);
%     x1{i}(end+1) = x2{idx}(1);
%     y1{i}(end+1) = y2{idx}(1);
% end
% 
% 
% for i = 1:(length(x2)-1)
%     val_aux = 0*[1:length(x1)];
%     for j = 1:length(x1)
%         aux = (cross([x2{i}(end) y2{i}(end) 0],[x1{j}(1) y1{j}(1) 0]));
%         if(aux(3)<0)
%             aux(3) = inf;
%         end
%         val_aux(j)  = aux(3);
%     end
%     [~, idx]=min(val_aux);
%     x2{i}(end+1) = x1{idx}(1);
%     y2{i}(end+1) = y1{idx}(1);
% end
% 
% 
% x = [];
% y = [];
% for i = 1:length(x1)
%     x = [x NaN x1{i}];
%     y = [y NaN y1{i}];
% end
% for i = 1:length(x2)
%     x = [x NaN x2{i}];
%     y = [y NaN y2{i}];
% end
% end
end
%% splitNaN
function [ out ] = splitNaN( a )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if ~isnan(a(end))
  a(end+1)=nan;
end
idx2=find(isnan(a));
idx1=[1 idx2(1:end-1)+1];
n=numel(idx1);
%out=cell(n,1);
j=1;
for k=1:n
    if(idx2(k)-1-idx1(k)>0)
      out{j}=a(idx1(k):idx2(k)-1);
      j=j+1;
    end
end

end

%% organize_points
function [x,y] = organize_points(xi,yi,N)
    if(~exist('N','var'))
        N= length(xi);
    end
    t = linspace(0,2*pi, N);
    val = xi*0;
    for i = 1:length(t)
        v = [sin(t(i)) cos(t(i))]';
        val = xi*0;
        for j=1:length(xi)
            val(j) = [xi(j) yi(j)]*v/norm([xi(j) yi(j)]);
        end
        [~,idx]=max(val);
        x(i) = xi(idx);
        y(i) = yi(idx);
        
    end
end