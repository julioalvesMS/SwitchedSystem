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

% Desired DC-DC conVsrter to use
% Options can be found in the system directory:
%   buck
%   boost
%   buck_boost
circuit = buck_boost(R, Ro, Co, L);

% Lambda used to create the mean system with wich the controle will be
% designed
% Value must be a Vsctor where:
%   sum(lambda) = 1
lambda = [0.5 0.5];

%% Prepare Data

run load_circuit_sys

%%

%Parâmetros
Ts = 5/10e4;
Tint = Ts/100;
Vo = 45;
ts = 0.25;
%
do = Vo/(Vo+Vs);
Po = Vo^2/Ro;
ILo = Po/Vs;
%Função de transferência do Boost
Kb = -ILo/Co;
z = (1/L)*(R-Vs/ILo);
a1 = R/L+1/(Co*Ro);
a2 = (1/(L*Co))*(R/Ro+(1-do)^2);

G = tf(Kb*[1 z],[1 a1 a2]);

bode(G)
%
%Para o Controlador:
Wcg0 = 600;%frequencia de cruzamento desejada
MF0 = 50*(pi/180);%margem de fase desejada
FaseG = angle(evalfr(G,Wcg0*1i));
Fo = MF0-FaseG-pi;  %passo 1
gama = -4.9*pi/180; %passo 2
Fm = Fo-gama;       %passo 3
% Projeto do Lead em casacata com o integrador
T1 = tan(gama+pi/2)/Wcg0; %passo 4
alfa = (1-sin(Fm))/(1+sin(Fm)) %passo 5
T2 = 1/(sqrt(alfa)*Wcg0); %passo 6
K0 = (Wcg0*sqrt(Wcg0^2+1/(alfa^2*T2^2)))/(sqrt(Wcg0^2+1/T2^2)*sqrt(Wcg0^2+1/T1^2)*abs(evalfr(G,Wcg0*1i)));  %passo 7
%controlador na forma continua:
K = tf(K0*conv([1 1/T1],[1 1/T2]),[1 1/(alfa*T2) 0]);

Kz = c2d(K,Ts,'tustin');
[nK,dK] = tfdata(Kz,'v');

figure
bode(K)
%
%Análise da planta
figure
margin(K*G)
%
figure
krl = 0:0.001:1;
rlocus(K*G,krl)
figure
pzmap(feedback(K*G,1))
%
%resposta ao degrau
figure
step(feedback(K*G,1))
%
%filtros
fc = 5000;
wc = 2*pi*fc;