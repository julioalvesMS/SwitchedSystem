%% Initial Setup
clear; clc; close all;

load('data/LC_01_22_grace.mat')

% folders to create
image_folder = 'images';
cache_folder = 'tmp/cache';

[~,~]=mkdir(image_folder);
image_folder = strcat(image_folder, '/');
[~,~]=mkdir(cache_folder);
cache_folder = strcat(cache_folder, '/');

addpath(genpath('functions'))
addpath(genpath('simulations'))
addpath(genpath('system_models'))
addpath(genpath('scripts'))


plot_compression_rate = 1e0;


Simulink.fileGenControl('set', 'CacheFolder', cache_folder);

%% System Specifications

run system_specifications

%% Sistema

% Desired Theorem to use
% Theorems defines
%   1 - Cost
%   2 - H2
%   3 - Hinf
opt_minimization = 2;

% Disturbances to be applied during simulations
% Options
%   disturbance_Vin_enable - Enable step disturbance in the input voltage
%   disturbance_Ro_enable - Enable step disturbance in the load resistance
disturbance_Vin_enable = false;
disturbance_Ro_enable = false;


% Desired DC-DC converter to use
% Options can be found  in the system directory:
%   buck
%   boost
%   buck_boost
%   buck_boost_non_inverting
circuit = buck(R, Ro, Co, L);

Vref = circuit.limit_cycle_voltage;

run load_circuit_sys

simulation_duration = 1;

%% Generate Limit Cycle

E = dsys.E;
H = dsys.H;
G = dsys.G;

Ed = zeros(2,2,dsys.N);
Hd = zeros(2,1,dsys.N);
Gd = zeros(2,1,dsys.N);
for i=length(dsys.N)
    Ed(:,:,i) = dsys.E{i};
    Gd(:,:,i) = dsys.G{i};
    Hd(:,:,i) = dsys.H{i};
end

xe = [0; Vref];
Gamma = circuit.limit_cycle_gamma;

[cand, kappa] = find_cycles(dsys, xe, Gamma);

cycle = find_limit_cycle(dsys, kappa, cand, opt_minimization);

%% Prepare Data

run circuit_disturbance


%% Simulate Converter 

% Get model to simulate
model = circuit.limit_cycle_simulink;

load_system(model);

% Convert the truct used to represent the space state to double, so it
% can be used in the simulink
[Ad, Bd, Cd, Dd, Qd, Ld] = gss2double(dsys);    % Discrete

ell = zeros(length(sys.A{1}), kappa, sys.N);
for i=1:dsys.N
    ell(:,:,i) = cell2mat({cycle.ell{i,:}});
end

P = zeros(length(sys.A{1}), length(sys.A{1}), kappa);
for i=1:kappa
    P(:,:,i) = cycle.lyap.P{i};
end

xe_h = cell2mat(cycle.xe_h);

rho = cycle.lyap.optval;

% Creates a bus, wich will be used in the simulink to simplify the
% model
SystemDataBus = create_bus_LimitCycleDataBus(Ad, Bd, P, Qd, sys.N, xe_h, ell, rho, Ed, Gd, Hd);

% Run simulation
sim(model, simulation_duration);

i = 1;
% Store only samples of the data, this will be made in order to save
% memory use
sim_out(i).IL = downsample_timeseries(logsout.get('IL').Values, plot_compression_rate);
sim_out(i).Vout = downsample_timeseries(logsout.get('Vout').Values, plot_compression_rate);


%% Analysis

plot_voltage_current(sim_out, circuit.name, image_folder);

plot_voltage_time(sim_out, circuit.name, image_folder);

plot_current_time(sim_out, circuit.name, image_folder);

if disturbance_Ro_enable == 1
    plot_disturbance_voltage_time(sim_out, disturbance_Ro_time, circuit.name, image_folder);
end

if disturbance_Vin_enable == 1
    plot_disturbance_voltage_time(sim_out, disturbance_Vin_time, circuit.name, image_folder);
end

%%

A = dsys.A;
Q = dsys.Q;
C = dsys.C;
D = dsys.D;
b = dsys.b;
x0 = dsys.x0;
N = length(dsys.A);
n = length(dsys.A{1});



xe_h = cycle.xe_h;
P = cycle.lyap.P;
ell = cycle.ell;

T = 1e-1;

switch opt_minimization
    case 1
        sigma = @(i_,x_) sigma_h2(A, ell, xe_h, P, Q, i_, x_(:,i_-1));
        wf = @(i_) 0;
    case 2
        sigma = @(i_,x_) sigma_h2(A, ell, xe_h, P, Q, i_, x_(:,i_-1));
        wf = @(i_) 0;
        m=cycle.lyap.m;
        x0 = A{m}*xe_h{kappa}+b{m}+H{m};
    case 3
        rho = cycle.lyap.optval;
        sigma = @(i_,x_) sigma_hinf(A, ell, xe_h, P, Q, H, G, rho, i_, x_(:,i_-1));
        wf = @(i_) 0.04*sin(2*pi*60*T*i_*(i_<.05/T));
        x0 = xe_h{1};
end


K=floor(simulation_duration/Ts); %[s]/T

x=zeros(n,K);
x(:,1) =x0;
xi=zeros(n,K);

for i=2:(K+1)
    s(i-1)= sigma(i,x);
    w(i-1) = wf(i);
    x(:,i) = A{s(i-1)}*x(:,i-1)+b{s(i-1)}+H{s(i-1)}*w(i-1) ;
    xi(:,i-1) = (x(:,i-1)-xe_h{mod(i-2,kappa)+1});
    ze(:,i-1) = C{s(i-1)}*xi(:,i-1)+D{s(i-1)}*w(i-1);
    v(i-1) = xi(:,i-1)'*P{mod(i-2,kappa)+1}*xi(:,i-1);
end

if(opt_minimization<2)   
        cost = trace(ze*ze');
elseif(opt_minimization<3)   
        cost = trace(ze*ze'+G{m}'*G{m});
else
        cost = trace(ze*ze')/trace(w*w')
end


f=figure;
f.Position(3:4) = [560 400];
hold on
plot(CO_Buck.IL, CO_Buck.Vof)
plot(x(1,:),x(2,:),'LineWidth',.2,'Marker','.','MarkerSize',10)
for i=1:kappa
    scatter(xe_h{i}(1),xe_h{i}(2),'LineWidth',7,'Marker','*')
end
xlabel('$x_1[n]$','interpreter','latex','FontSize',12)
ylabel('$x_2[n]$','interpreter','latex','FontName','Times-Roman','FontSize',12)
%set(gca,'FontSize',15)


k= 0:(length(s)-1);
f=figure;
f.Position(3:4) = [560 240];
plot(k, xi','LineWidth',.2,'Marker','.','MarkerSize',10)
xlabel('$n$','interpreter','latex','FontName','Times-Roman','FontSize',12)
ylabel('$\xi[n]$','interpreter','latex','FontName','Times-Roman','FontSize',12)
ylim([-1 15])
xlim([0 K])
%set(gca,'FontSize',15)

f=figure;
f.Position(3:4) = [560 240];
stairs(k,s,'LineWidth',.2,'Marker','.','MarkerSize',10)
yticks([1 2])
ylim([.8 2.2])
xlim([0 K])
xlabel('$n$','interpreter','latex','FontName','Times-Roman','FontSize',12)
ylabel('$\sigma[n]$','interpreter','latex','FontName','Times-Roman','FontSize',12)
%axis([0 K 0.95 2.05])
%set(gca,'FontSize',15)

%end


function [s] = sigma_h2(A,ell,xe,P,Q,n,x) 
    kappa= length(P);
    k =mod(n-2,kappa)+1;
    km=mod(n-1,kappa)+1;
    for i=1:length(A)
        val(i) = [x-xe{k};1]'*([A{i}'*P{km}*A{i}-P{k}+Q{i}  A{i}'*P{km}*ell{i,k};  ell{i,k}'*P{km}*A{i} ell{i,k}'*P{km}*ell{i,k}]*[x-xe{k};1]);
    end
    [~,s] = min(val);
end

function [s] = sigma_hinf(A,ell,xe,P,E,H,G,rho,n,x) 
    kappa= length(P);
    k =mod(n-2,kappa)+1;
    km=mod(n-1,kappa)+1;
    for i=1:length(A)
        Lcal = [A{i}'*P{km}*A{i}-P{k}      A{i}'*P{km}*ell{i,k};
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

