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


plot_compression_rate = 1e3;


Simulink.fileGenControl('set', 'CacheFolder', cache_folder);

%% System Specifications

run system_specifications

%% Simulation Parametersr

% Desired DC-DC converter to use
% Options can be found in the system directory:
%   buck
%   boost
%   buck_boost
circuit = buck_boost(R, Ro, Co, L);

%% Prepare Data

run load_circuit_sys


% Lambda used to create the mean system with wich the controle will be
% designed
% Value must be a vector where:
%   sum(lambda) = 1
lambda = generate_lambda_voltage(sys, 100);
[Alamb, Blamb, ~] = calc_sys_lambda(sys, lambda);
Blamb = Blamb*Vs;

C = [0 1];
D = 0;

sys_ss = ss(Alamb, Blamb, C, D);

Cs = 3.28/s + 0.00919;

Gss = tf(sys_ss);

%Parâmetros
Vo = 100;
d = Vo/(Vo+Vs);
G1 = Vs*((d*L/R)*s - (1-d)^2)/((1-d)^2*(L*Co*s^2 + L*s/R + (1-d)^2))


%%

%Parâmetros
Vo = 100;
%
%syms L Ts Ro Co R Vs Vo s
do = Vo/(Vo+Vs);
Po = Vo^2/Ro;
ILo = Po/Vs;
VdV = Vs - V;
%Função de transferência do Boost
Kb = -Vgv/do^2;
z = -(L*ILo)/VgV;
a1 = (L*Co)/do^2;
a2 = (1/(L*Co))*(R/Ro+(1-do)^2);

G2_t = Kb*(s*z + 1)/(a2*s^2 + a1*s + 1)


%%

%Parâmetros
Vo = 100;
%
%syms L Ts Ro Co R Vs Vo s
do = Vo/(Vo+Vs);
Po = Vo^2/Ro;
ILo = Po/Vs;
%Função de transferência do Boost
Kb = -ILo/Co;
z = (1/L)*(R-Vs/ILo);
a1 = R/L+1/(Co*Ro);
a2 = (1/(L*Co))*(R/Ro+(1-do)^2);

G2 = Kb*(s + z)/(s^2 + a1*s +a2)

%G2 = tf(Kb*[1 z],[1 a1 a2]);

%%
close all
G = G2;

% Wcg0 = 3.5e2;%frequencia de cruzamento desejada
% MF0 = 65*(pi/180);%margem de fase desejada
% Wcg0 = 10.9;%frequencia de cruzamento desejada
% MF0 = 108*(pi/180);%margem de fase desejada
Wcg0 = 190;%frequencia de cruzamento desejada
MF0 = 60*(pi/180);%margem de fase desejada
FaseG = angle(evalfr(G,Wcg0*1i));
Fm = MF0-FaseG-pi/2;  %passo 1
Ti = tan(Fm)/Wcg0;    %passo 2
Kp = 1/(sqrt(1+(Ti*Wcg0)^-2)*abs(evalfr(G,Wcg0*1i))); %passo 3
K = tf(Kp*[1 1/Ti],[1 0]);
F = feedback(K*G,1);

figure
margin(K*G)
pid(K)

return
%%

kp = circuit.pwm_pid_kp;
ki = circuit.pwm_pid_ki;

Kr = kp + ki/s;

Fr = feedback(Kr*G,1);

figure
margin(Kr*G)