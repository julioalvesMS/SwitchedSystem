%% Initial Setup
clear; clc; close all;

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

Simulink.fileGenControl('set', 'CacheFolder', cache_folder);

%% System Specifications

run system_specifications

%% Simulation Parametersr

% Desired DC-DC converter to use
% Options can be found in the system directory:
%   buck
%   boost
%   buck_boost
%   buck_boost_non_inverting
circuit = buck_boost_non_inverting(R, Ro, Co, L);

%% Prepare Data

run load_circuit_sys

%% System Equilibrium Points

step = 0.001;
sequence = 0:step:1;

test_lambdas = zeros(498842, 3);
index = 1;
for lambda1=sequence
    for lambda2=0:step:(1-lambda1)
        lambda3 = 1 - lambda2 - lambda1;
        test_lambdas(index,1) = lambda1;
        test_lambdas(index,2) = lambda2;
        test_lambdas(index,3) = lambda3;
        index = index+1;
    end
end

% Calculate the equilibrium point for each lambda sample 
for i=size(test_lambdas, 1):-1:1
    [Al, Bl, ~] = calc_sys_lambda(sys, test_lambdas(i,:));
    xe = -Al\Bl*sys.U;

    equilibrium(i,:) = xe;
end

%% Analysis

figure
plot(equilibrium(:,1), equilibrium(:,2), 'd')
ylabel('v_o [V]');
xlabel('IL [A]');
saveas(gcf, strcat(image_folder, 'equilibrium'), 'eps')

% plot_voltage_current_equilibrium(equilibrium, circuit.name, image_folder);
% plot_voltage_lambda_equilibrium(equilibrium, sample_lambdas, circuit.name, image_folder);

%%
X = equilibrium(:,2);
Y = equilibrium(:,1);

X0 = min(X);
[XF, iXF] = max(X);
YXF = Y(iXF);

Xstep = (XF - X0)/1e3;
Icy = find(Y < YXF);

Yc = Y(Icy);
Xc = X(Icy);

Ig = [];
Wg = [];
for x = X0:Xstep:(XF-Xstep)
    I = find(Xc > x & Xc < (x+Xstep));
    [~,iY] = min(Yc(I));
    
    Ig(end+1) = I(iY);
    Wg(end+1) = length(I);
end

Xmin = Xc(Ig);
Ymin = Yc(Ig);

figure
hold on
plot(X, Y, 'd')
plot(Xmin, Ymin, 'd')
hold off

cftool(Xmin, Ymin, [], Wg)

%%

XP = X0:(XF-X0)/1e4:XF;

f=fit(Xmin,Ymin,'poly8', 'Normalize', 'on', 'Weight', Wg);

figure
hold on
plot(equilibrium(:,2), equilibrium(:,1), 'd')
plot(XP, f(XP), '--r')
hold off


figure
hold on
plot(equilibrium(:,1), equilibrium(:,2), 'd')
plot(f(XP), XP, '--r')
hold off