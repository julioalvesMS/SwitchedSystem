%% Initial Setup
clear; clc; close all;

% folders to create
root_image_folder = 'images';
cache_folder = 'tmp/cache';

[~,~]=mkdir(root_image_folder);
root_image_folder = strcat(root_image_folder, '/');
[~,~]=mkdir(cache_folder);
cache_folder = strcat(cache_folder, '/');

addpath(genpath('functions'))
addpath(genpath('simulations'))  
addpath(genpath('system_models'))
addpath(genpath('scripts'))
addpath(genpath('data'))


plot_compression_rate = 1e0;


Simulink.fileGenControl('set', 'CacheFolder', cache_folder);


image_folder = strcat(root_image_folder, '/Article');
[~,~]=mkdir(image_folder);
image_folder = strcat(image_folder, '/');

image_format = 'epsc';

continuous_color = [0, 0.4470, 0.7410];
discrete_color = [0.8500, 0.3250, 0.0980];
classic_color = [0.9290, 0.6940, 0.1250];

root_data_folder = 'data';
root_data_mat_folder = 'data_mat';
data_folder = "CURVAS_2020_08_25";

%% System Specifications

run system_specifications

run circuit_disturbance

% Desired DC-DC converter to use
% Options can be found in the system directory:
%   buck
%   boost
%   buck_boost
%   buck_boost_non_inverting
circuit = buck_boost(R, Ro, Co, L);

run load_circuit_sys

%% End of simulations

file = fullfile(root_data_mat_folder, 'all_simulations_article.mat');

load(file)

%% Load experiment data

folder_path = fullfile(root_data_folder, data_folder, '*.dat');

data_files = dir(folder_path);
N = length(data_files);
 
for i = 1:N
    file = data_files(i);
    path = fullfile(file.folder, file.name);
    load(path);
end

ratio = 1;

% Degrau - Clássico - 100V
data = {};
data_C1 = C1STEP_PWM_100V00001;
data_C4 = C4STEP_PWM_100V00001;
data.t = data_C1(1:ratio:end,1);
data.Vout = timeseries(data_C1(1:ratio:end,2), data.t);
data.Vref = timeseries(100*ones(size(data.Vout)), data.t);
data.IL = timeseries(data_C4(1:ratio:end,2), data.t);
exp_step_pwm_100 = data;

% Degrau - Contínuo 1 - 100V
data = {};
data_C1 = C1STEP_CON1_100V00001;
data_C4 = C4STEP_CON1_100V00001;
data.t = data_C1(1:ratio:end,1);
data.Vout = timeseries(data_C1(1:ratio:end,2), data.t);
data.Vref = timeseries(100*ones(size(data.Vout)), data.t);
data.IL = timeseries(data_C4(1:ratio:end,2), data.t);
exp_step_con1_100 = data;

% Degrau - Contínuo 2 - 100V
data = {};
data_C1 = C1STEP_CON2_100V00001;
data_C4 = C4STEP_CON2_100V00001;
data.t = data_C1(1:ratio:end,1);
data.Vout = timeseries(data_C1(1:ratio:end,2), data.t);
data.Vref = timeseries(100*ones(size(data.Vout)), data.t);
data.IL = timeseries(data_C4(1:ratio:end,2), data.t);
exp_step_con2_100 = data;

% Degrau - Discreto - 100V
data = {};
data_C1 = C1STEP_DIS_100V00101;
data_C4 = C4STEP_DIS_100V00101;
data.t = data_C1(1:ratio:end,1);
data.Vout = timeseries(data_C1(1:ratio:end,2), data.t);
data.Vref = timeseries(100*ones(size(data.Vout)), data.t);
data.IL = timeseries(data_C4(1:ratio:end,2), data.t);
exp_step_dis_100 = data;




% Degrau Carga - Clássico - 100V
data = {};
data_C1 = C1LOAD_PWM_100V00001;
data_C4 = C4LOAD_PWM_100V00001;
data.t = data_C1(1:ratio:end,1);
data.Vout = timeseries(data_C1(1:ratio:end,2), data.t);
data.Vref = timeseries(100*ones(size(data.Vout)), data.t);
data.IL = timeseries(data_C4(1:ratio:end,2), data.t);
exp_load_pwm_100 = data;

% Degrau Carga - Contínuo 1 - sem Correção de corente - 100V
data = {};
data_C1 = C1LOAD_CON1_noC_100V00000;
data_C4 = C4LOAD_CON1_noC_100V00000;
data.t = data_C1(1:ratio:end,1);
data.Vout = timeseries(data_C1(1:ratio:end,2), data.t);
data.Vref = timeseries(100*ones(size(data.Vout)), data.t);
data.IL = timeseries(data_C4(1:ratio:end,2), data.t);
exp_load_con1_noC_100 = data;

% Degrau Carga - Contínuo 1 - com Correção de corente - 100V
data = {};
data_C1 = C1LOAD_CON1_C_100V00001;
data_C4 = C4LOAD_CON1_C_100V00001;
data.t = data_C1(1:ratio:end,1);
data.Vout = timeseries(data_C1(1:ratio:end,2), data.t);
data.Vref = timeseries(100*ones(size(data.Vout)), data.t);
data.IL = timeseries(data_C4(1:ratio:end,2), data.t);
exp_load_con1_C_100 = data;

% Degrau Carga - Contínuo 2 - sem Correção de corente - 100V
data = {};
data_C1 = C1LOAD_CON2_noC_100V00000;
data_C4 = C4LOAD_CON2_noC_100V00000;
data.t = data_C1(1:ratio:end,1);
data.Vout = timeseries(data_C1(1:ratio:end,2), data.t);
data.Vref = timeseries(100*ones(size(data.Vout)), data.t);
data.IL = timeseries(data_C4(1:ratio:end,2), data.t);
exp_load_con2_noC_100 = data;

% Degrau Carga - Contínuo 2 - com Correção de corente - 100V
data = {};
data_C1 = C1LOAD_CON2_C_100V00001;
data_C4 = C4LOAD_CON2_C_100V00001;
data.t = data_C1(1:ratio:end,1);
data.Vout = timeseries(data_C1(1:ratio:end,2), data.t);
data.Vref = timeseries(100*ones(size(data.Vout)), data.t);
data.IL = timeseries(data_C4(1:ratio:end,2), data.t);
exp_load_con2_C_100 = data;

% Degrau Carga - Discreto - sem Correção de corente - 100V
data = {};
data_C1 = C1LOAD_DIS_noC_100V00000;
data_C4 = C4LOAD_DIS_noC_100V00000;
data.t = data_C1(1:ratio:end,1);
data.Vout = timeseries(data_C1(1:ratio:end,2), data.t);
data.Vref = timeseries(100*ones(size(data.Vout)), data.t);
data.IL = timeseries(data_C4(1:ratio:end,2), data.t);
exp_load_dis_noC_100 = data;

% Degrau Carga - Discreto - com Correção de corente - 100V
data = {};
data_C1 = C1LOAD_DIS_C_100V00001;
data_C4 = C4LOAD_DIS_C_100V00001;
data.t = data_C1(1:ratio:end,1);
data.Vout = timeseries(data_C1(1:ratio:end,2), data.t);
data.Vref = timeseries(100*ones(size(data.Vout)), data.t);
data.IL = timeseries(data_C4(1:ratio:end,2), data.t);
exp_load_dis_C_100 = data;



for i = 1:N
    var = split(data_files(i).name, '.');
    clear(var{1})
end

% Remove noise
exp_load_con1_C_100.Vout.Data = smooth(exp_load_con1_C_100.Vout.Time,exp_load_con1_C_100.Vout.Data,6,'rloess');
exp_load_con2_C_100.Vout.Data = smooth(exp_load_con2_C_100.Vout.Time,exp_load_con2_C_100.Vout.Data,6,'rloess');
exp_load_dis_C_100.Vout.Data = smooth(exp_load_dis_C_100.Vout.Time,exp_load_dis_C_100.Vout.Data,6,'rloess');
exp_load_pwm_100.Vout.Data = smooth(exp_load_pwm_100.Vout.Time,exp_load_pwm_100.Vout.Data,6,'rloess');

exp_step_con1_100.Vout.Data = smooth(exp_step_con1_100.Vout.Time,exp_step_con1_100.Vout.Data,6,'rloess');
exp_step_con2_100.Vout.Data = smooth(exp_step_con2_100.Vout.Time,exp_step_con2_100.Vout.Data,6,'rloess');
exp_step_dis_100.Vout.Data = smooth(exp_step_dis_100.Vout.Time,exp_step_dis_100.Vout.Data,6,'rloess');
exp_step_pwm_100.Vout.Data = smooth(exp_step_pwm_100.Vout.Time,exp_step_pwm_100.Vout.Data,6,'rloess');


%% CSV
st = 2;

con1 = csvread('exp_frequency_con1_noMH_noC.csv');
con2 = csvread('exp_frequency_con2_noMH_noC.csv');
dis = csvread('exp_frequency_disc_noMH_noC.csv');
data.Vref = dis(st:end,1);
data.Fcon1 = con1(st:end,2);
data.Fcon2 = con2(st:end,2);
data.Fdis = dis(st:end,2);
data.Vcon1 = con1(st:end,3);
data.Vcon2 = con2(st:end,3);
data.Vdis = dis(st:end,3);
exp_frequency_noMH_noC = data;

con1 = csvread('exp_frequency_con1_noMH_C.csv');
con2 = csvread('exp_frequency_con2_noMH_C.csv');
dis = csvread('exp_frequency_disc_noMH_C.csv');
data.Vref = dis(st:end,1);
data.Fcon1 = con1(st:end,2);
data.Fcon2 = con2(st:end,2);
data.Fdis = dis(st:end,2);
data.Vcon1 = con1(st:end,3);
data.Vcon2 = con2(st:end,3);
data.Vdis = dis(st:end,3);
exp_frequency_noMH_C = data;

con1 = csvread('exp_frequency_con1_MH_noC.csv');
con2 = csvread('exp_frequency_con2_MH_noC.csv');
dis = csvread('exp_frequency_disc_MH_noC.csv');
data.Vref = dis(st:end,1);
data.Fcon1 = con1(st:end,2);
data.Fcon2 = con2(st:end,2);
data.Fdis = dis(st:end,2);
data.Vcon1 = con1(st:end,3);
data.Vcon2 = con2(st:end,3);
data.Vdis = dis(st:end,3);
exp_frequency_MH_noC = data;

con1 = csvread('exp_frequency_con1_MH_C.csv');
con2 = csvread('exp_frequency_con2_MH_C.csv');
dis = csvread('exp_frequency_disc_MH_C.csv');
data.Vref = dis(st:end,1);
data.Fcon1 = con1(st:end,2);
data.Fcon2 = con2(st:end,2);
data.Fdis = dis(st:end,2);
data.Vcon1 = con1(st:end,3);
data.Vcon2 = con2(st:end,3);
% data.Vdis = dis(st:end,3);
exp_frequency_MH_C = data;

%% Figure generation


%% Equilibrium Region

[sample_lambdas, equilibrium] = generate_sample_points(sys);


figure
plot(equilibrium(:,1), equilibrium(:,2), 'black','LineWidth',1.5)
ylabel('Output voltage (V)');
xlabel('i_L (A)');
ax = gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;
set(gcf,'renderer','Painters')
saveas(gcf, strcat(image_folder, 'equilibrium'), image_format)


%% Frequency - Frequency Variation at 1MHz - Ideal Switches

Vref = continuous_1_frequency_ideal_1M.sim_out.Vref;
Freq = continuous_1_frequency_ideal_1M.sim_out.F;
Tref = continuous_1_frequency_ideal_1M.variable_reference_step_period;
[continuous_1_F, continuous_1_FV] = map_voltage_frequency(Vref, Freq, Tref);

Vref = continuous_2_frequency_ideal_1M.sim_out.Vref;
Freq = continuous_2_frequency_ideal_1M.sim_out.F;
Tref = continuous_2_frequency_ideal_1M.variable_reference_step_period;
[continuous_2_F, continuous_2_FV] = map_voltage_frequency(Vref, Freq, Tref);

Vref = discrete_frequency_ideal_1M.sim_out.Vref;
Freq = discrete_frequency_ideal_1M.sim_out.F;
Tref = discrete_frequency_ideal_1M.variable_reference_step_period;
[discrete_F, discrete_FV] = map_voltage_frequency(Vref, Freq, Tref);

figure
hold all
plot(continuous_2_FV, continuous_2_F/1e3, 'DisplayName', 'QNS Controller', 'LineWidth',1.5)
plot(continuous_1_FV, continuous_1_F/1e3, 'DisplayName', 'RNS Controller', 'LineWidth',1.5)
plot(discrete_FV, discrete_F/1e3, 'DisplayName', 'RS Controller', 'LineWidth',1.5)
hold off
ylabel('Switching frequency (kHz)')
xlabel('Output voltage (V)')
le = legend;
le.FontSize = 15;
set(le, 'Position',[0.514880959618659 0.226984130390106 0.333928564190864 0.165476185934884]);
ax = gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;
set(gcf,'renderer','Painters')
saveas(gcf, strcat(image_folder, 'frequency_variation_1MHz'), image_format);


%% Frequency - Frequency Variation at 200kHz - Ideal Switches

Vref = continuous_1_frequency_ideal_200.sim_out.Vref;
Freq = continuous_1_frequency_ideal_200.sim_out.F;
Tref = continuous_1_frequency_ideal_200.variable_reference_step_period;
[continuous_1_F, continuous_1_FV] = map_voltage_frequency(Vref, Freq, Tref);

Vref = continuous_2_frequency_ideal_200.sim_out.Vref;
Freq = continuous_2_frequency_ideal_200.sim_out.F;
Tref = continuous_2_frequency_ideal_200.variable_reference_step_period;
[continuous_2_F, continuous_2_FV] = map_voltage_frequency(Vref, Freq, Tref);

Vref = discrete_frequency_ideal_200.sim_out.Vref;
Freq = discrete_frequency_ideal_200.sim_out.F;
Tref = discrete_frequency_ideal_200.variable_reference_step_period;
[discrete_F, discrete_FV] = map_voltage_frequency(Vref, Freq, Tref);

figure
hold all
plot(continuous_2_FV, continuous_2_F/1e3, 'DisplayName', 'QNS Controller', 'LineWidth',1.5)
plot(continuous_1_FV, continuous_1_F/1e3, 'DisplayName', 'RNS Controller', 'LineWidth',1.5)
plot(discrete_FV, discrete_F/1e3, 'DisplayName', 'RS Controller', 'LineWidth',1.5)
hold off
ylabel('Switching frequency (kHz)')
xlabel('Output voltage (V)')
le = legend;
le.FontSize = 15;
set(le, 'Position',[0.514880959618659 0.226984130390106 0.333928564190864 0.165476185934884]);
ax = gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;
set(gcf,'renderer','Painters')
saveas(gcf, strcat(image_folder, 'frequency_variation_200kHz'), image_format);


%% Frequency - Frequency Variation at 40kHz - Ideal Switches

Vref = continuous_1_frequency_ideal_40.sim_out.Vref;
Freq = continuous_1_frequency_ideal_40.sim_out.F;
Tref = continuous_1_frequency_ideal_40.variable_reference_step_period;
[continuous_1_F, continuous_1_FV] = map_voltage_frequency(Vref, Freq, Tref);

Vref = continuous_2_frequency_ideal_40.sim_out.Vref;
Freq = continuous_2_frequency_ideal_40.sim_out.F;
Tref = continuous_2_frequency_ideal_40.variable_reference_step_period;
[continuous_2_F, continuous_2_FV] = map_voltage_frequency(Vref, Freq, Tref);

Vref = discrete_frequency_ideal_40.sim_out.Vref;
Freq = discrete_frequency_ideal_40.sim_out.F;
Tref = discrete_frequency_ideal_40.variable_reference_step_period;
[discrete_F, discrete_FV] = map_voltage_frequency(Vref, Freq, Tref);

figure
hold all
plot(continuous_2_FV, continuous_2_F/1e3, 'DisplayName', 'QNS Controller', 'LineWidth',1.5)
plot(continuous_1_FV, continuous_1_F/1e3, 'DisplayName', 'RNS Controller', 'LineWidth',1.5)
plot(discrete_FV, discrete_F/1e3, 'DisplayName', 'RS Controller', 'LineWidth',1.5)
hold off
ylabel('Switching frequency (kHz)')
xlabel('Output voltage (V)')
le = legend;
le.FontSize = 15;
set(le, 'Position',[0.514880959618659 0.226984130390106 0.333928564190864 0.165476185934884]);
ax = gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;
set(gcf,'renderer','Painters')
saveas(gcf, strcat(image_folder, 'frequency_variation_40kHz'), image_format);


%% Voltage Error - Frequency Variation at 1MHz - Ideal Switches

Vref = continuous_1_frequency_ideal_1M.sim_out.Vref;
Verr = continuous_1_frequency_ideal_1M.sim_out.Verr;
Tref = continuous_1_frequency_ideal_1M.variable_reference_step_period;
[continuous_1_Er, continuous_1_Vr] = map_voltage_frequency(Vref, Verr, Tref);

Vref = continuous_2_frequency_ideal_1M.sim_out.Vref;
Verr = continuous_2_frequency_ideal_1M.sim_out.Verr;
Tref = continuous_2_frequency_ideal_1M.variable_reference_step_period;
[continuous_2_Er, continuous_2_Vr] = map_voltage_frequency(Vref, Verr, Tref);

Vref = discrete_frequency_ideal_1M.sim_out.Vref;
Verr = discrete_frequency_ideal_1M.sim_out.Verr;
Tref = discrete_frequency_ideal_1M.variable_reference_step_period;
[discrete_Er, discrete_Vr] = map_voltage_frequency(Vref, Verr, Tref);

figure
hold all
plot(continuous_2_Vr, continuous_2_Er, 'DisplayName', 'QNS Controller', 'LineWidth',1.5)
plot(continuous_1_Vr, continuous_1_Er, 'DisplayName', 'RNS Controller', 'LineWidth',1.5)
plot(discrete_Vr, discrete_Er, 'DisplayName', 'RS Controller', 'LineWidth',1.5)
hold off
grid
ylabel('Voltage Error (V)')
xlabel('Output voltage (V)')
le = legend;
le.FontSize = 15;
% legend1 = legend;
% set(legend1, 'Position',[0.514880959618659 0.176984130390106 0.333928564190864 0.165476185934884]);
%set(gcf, 'Position',  [100, 100, 1600, 780])
ax = gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;
set(gcf,'renderer','Painters')
saveas(gcf, strcat(image_folder, 'voltage_error_1MHz'), image_format);


%% Voltage Error - Frequency Variation at 200kHz - Ideal Switches

Vref = continuous_1_frequency_ideal_200.sim_out.Vref;
Verr = continuous_1_frequency_ideal_200.sim_out.Verr;
Tref = continuous_1_frequency_ideal_200.variable_reference_step_period;
[continuous_1_Er, continuous_1_Vr] = map_voltage_frequency(Vref, Verr, Tref);

Vref = continuous_2_frequency_ideal_200.sim_out.Vref;
Verr = continuous_2_frequency_ideal_200.sim_out.Verr;
Tref = continuous_2_frequency_ideal_200.variable_reference_step_period;
[continuous_2_Er, continuous_2_Vr] = map_voltage_frequency(Vref, Verr, Tref);

Vref = discrete_frequency_ideal_200.sim_out.Vref;
Verr = discrete_frequency_ideal_200.sim_out.Verr;
Tref = discrete_frequency_ideal_200.variable_reference_step_period;
[discrete_Er, discrete_Vr] = map_voltage_frequency(Vref, Verr, Tref);

figure
hold all
plot(continuous_2_Vr, continuous_2_Er, 'DisplayName', 'QNS Controller', 'LineWidth',1.5)
plot(continuous_1_Vr, continuous_1_Er, 'DisplayName', 'RNS Controller', 'LineWidth',1.5)
plot(discrete_Vr, discrete_Er, 'DisplayName', 'RS Controller', 'LineWidth',1.5)
hold off
grid
ylabel('Voltage Error (V)')
xlabel('Output voltage (V)')
le = legend;
le.FontSize = 15;
% legend1 = legend;
% set(legend1, 'Position',[0.514880959618659 0.176984130390106 0.333928564190864 0.165476185934884]);
%set(gcf, 'Position',  [100, 100, 1600, 780])
ax = gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;
set(gcf,'renderer','Painters')
saveas(gcf, strcat(image_folder, 'voltage_error_200kHz'), image_format);


%% Voltage Error - Frequency Variation at 40kHz - Ideal Switches

Vref = continuous_1_frequency_ideal_40.sim_out.Vref;
Verr = continuous_1_frequency_ideal_40.sim_out.Verr;
Tref = continuous_1_frequency_ideal_40.variable_reference_step_period;
[continuous_1_Er, continuous_1_Vr] = map_voltage_frequency(Vref, Verr, Tref);

Vref = continuous_2_frequency_ideal_40.sim_out.Vref;
Verr = continuous_2_frequency_ideal_40.sim_out.Verr;
Tref = continuous_2_frequency_ideal_40.variable_reference_step_period;
[continuous_2_Er, continuous_2_Vr] = map_voltage_frequency(Vref, Verr, Tref);

Vref = discrete_frequency_ideal_40.sim_out.Vref;
Verr = discrete_frequency_ideal_40.sim_out.Verr;
Tref = discrete_frequency_ideal_40.variable_reference_step_period;
[discrete_Er, discrete_Vr] = map_voltage_frequency(Vref, Verr, Tref);

figure
hold all
plot(continuous_2_Vr, continuous_2_Er*100./continuous_2_Vr, 'DisplayName', 'QNS Controller', 'LineWidth',1.5)
plot(continuous_1_Vr, continuous_1_Er*100./continuous_2_Vr, 'DisplayName', 'RNS Controller', 'LineWidth',1.5)
plot(discrete_Vr, discrete_Er*100./continuous_2_Vr, 'DisplayName', 'RS Controller', 'LineWidth',1.5)
hold off
grid
ylabel('Voltage Error (V)')
xlabel('Output voltage (V)')
legend
% legend1 = legend;
% set(legend1, 'Position',[0.514880959618659 0.176984130390106 0.333928564190864 0.165476185934884]);
%set(gcf, 'Position',  [100, 100, 1600, 780])
ax = gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;
set(gcf,'renderer','Painters')
saveas(gcf, strcat(image_folder, 'voltage_error_40kHz'), image_format);


%% Frequency - Real 40kHz - PI - Continuous Controller 1

Vref = continuous_1_frequency_ideal_40.sim_out.Vref;
Freq = continuous_1_frequency_ideal_40.sim_out.F;
Tref = continuous_1_frequency_ideal_40.variable_reference_step_period;
[switching_ideal_F, switching_ideal_FV] = map_voltage_frequency(Vref, Freq, Tref);

Vref = continuous_1_frequency_real_PI.sim_out.Vref;
Freq = continuous_1_frequency_real_PI.sim_out.F;
Tref = continuous_1_frequency_real_PI.variable_reference_step_period;
[switching_real_F, switching_real_FV] = map_voltage_frequency(Vref, Freq, Tref);

switching_classic_FV = switching_real_FV;
switching_classic_F = 40e3*ones(size((switching_classic_FV)));

exp_FV = exp_frequency_noMH_C.Vref;
exp_F = exp_frequency_noMH_C.Fcon1*2;

figure
hold all
plot(switching_ideal_FV, switching_ideal_F/1e3, '--', 'DisplayName', 'Ideal Switches', 'LineWidth',1.5,'color','blue')
plot(switching_real_FV, switching_real_F/1e3, '-x', 'DisplayName', 'Real Gates', 'LineWidth',1.5,'color','blue')
plot(switching_classic_FV, switching_classic_F/1e3, '--', 'DisplayName', 'Equivalent PWM', 'LineWidth',1.5,'color','black')
plot(exp_FV, exp_F, '-o', 'DisplayName', 'Experimental Result', 'LineWidth',1.5,'color', 'red')
hold off
le = legend;
le.FontSize = 15;
set(le, 'Position',[0.514880959618659 0.226984130390106 0.333928564190864 0.165476185934884]);
ylabel('Switching frequency (kHz)')
xlabel('Output voltage (V)')
%set(gcf, 'Position',  [100, 100, 1600, 780])
ax = gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;
set(gcf,'renderer','Painters')
saveas(gcf, strcat(image_folder, 'frequency_simulation_continuous_1'), image_format);


%% Frequency - Real 40kHz - PI - Continuous Controller 2

Vref = continuous_2_frequency_ideal_PI_40.sim_out.Vref;
Freq = continuous_2_frequency_ideal_PI_40.sim_out.F;
Tref = continuous_2_frequency_ideal_PI_40.variable_reference_step_period;
[switching_ideal_F, switching_ideal_FV] = map_voltage_frequency(Vref, Freq, Tref);

Vref = continuous_2_frequency_real_PI.sim_out.Vref;
Freq = continuous_2_frequency_real_PI.sim_out.F;
Tref = continuous_2_frequency_real_PI.variable_reference_step_period;
[switching_real_F, switching_real_FV] = map_voltage_frequency(Vref, Freq, Tref);

switching_classic_FV = switching_real_FV;
switching_classic_F = 40e3*ones(size((switching_classic_FV)));

exp_FV = exp_frequency_noMH_C.Vref;
exp_F = exp_frequency_noMH_C.Fcon2*2;

figure
hold all
plot(switching_ideal_FV, switching_ideal_F/1e3, '--', 'DisplayName', 'Ideal Switches', 'LineWidth',1.5,'color','blue')
plot(switching_real_FV, switching_real_F/1e3, '-x', 'DisplayName', 'Real Gates', 'LineWidth',1.5,'color','blue')
plot(switching_classic_FV, switching_classic_F/1e3, '--', 'DisplayName', 'Equivalent PWM', 'LineWidth',1.5,'color','black')
plot(exp_FV, exp_F, '-o', 'DisplayName', 'Experimental Result', 'LineWidth',1.5,'color', 'red')
hold off
le = legend;
le.FontSize = 15;
set(le, 'Position',[0.514880959618659 0.226984130390106 0.333928564190864 0.165476185934884]);
ylabel('Switching frequency (kHz)')
xlabel('Output voltage (V)')
%set(gcf, 'Position',  [100, 100, 1600, 780])
ax = gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;
set(gcf,'renderer','Painters')
saveas(gcf, strcat(image_folder, 'frequency_simulation_continuous_2'), image_format);


%% Frequency - Real 40kHz - PI - Discrete Controller

Vref = discrete_frequency_ideal_PI_40.sim_out.Vref;
Freq = discrete_frequency_ideal_PI_40.sim_out.F;
Tref = discrete_frequency_ideal_PI_40.variable_reference_step_period;
[switching_ideal_F, switching_ideal_FV] = map_voltage_frequency(Vref, Freq, Tref);

Vref = discrete_frequency_real_PI.sim_out.Vref;
Freq = discrete_frequency_real_PI.sim_out.F;
Tref = discrete_frequency_real_PI.variable_reference_step_period;
[switching_real_F, switching_real_FV] = map_voltage_frequency(Vref, Freq, Tref);

switching_classic_FV = switching_real_FV;
switching_classic_F = 40e3*ones(size((switching_classic_FV)));

exp_FV = exp_frequency_noMH_C.Vref;
exp_F = exp_frequency_noMH_C.Fdis*2;

figure
hold all
plot(switching_ideal_FV, switching_ideal_F/1e3, '--', 'DisplayName', 'Ideal Switches', 'LineWidth',1.5,'color','blue')
plot(switching_real_FV, switching_real_F/1e3, '-x', 'DisplayName', 'Real Gates', 'LineWidth',1.5,'color','blue')
plot(switching_classic_FV, switching_classic_F/1e3, '--', 'DisplayName', 'Equivalent PWM', 'LineWidth',1.5,'color','black')
plot(exp_FV, exp_F, '-o', 'DisplayName', 'Experimental Result', 'LineWidth',1.5,'color', 'red')
hold off
le = legend;
le.FontSize = 15;
set(le, 'Position',[0.514880959618659 0.226984130390106 0.333928564190864 0.165476185934884]);
ylabel('Switching frequency (kHz)')
xlabel('Output voltage (V)')
%set(gcf, 'Position',  [100, 100, 1600, 780])
ax = gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;
set(gcf,'renderer','Painters')
saveas(gcf, strcat(image_folder, 'frequency_simulation_discrete'), image_format);


%% Ripple

Vref = continuous_1_frequency_real_PI.sim_out.Vref;
Vrip = continuous_1_frequency_real_PI.sim_out.Ripple;
Tref = continuous_1_frequency_real_PI.variable_reference_step_period;
[continuous_1_R,continuous_1_RV] = map_voltage_frequency(Vref, Vrip, Tref);

Vref = continuous_2_frequency_real_PI.sim_out.Vref;
Vrip = continuous_2_frequency_real_PI.sim_out.Ripple;
Tref = continuous_2_frequency_real_PI.variable_reference_step_period;
[continuous_R,continuous_RV] = map_voltage_frequency(Vref, Vrip, Tref);

Vref = discrete_frequency_real_PI.sim_out.Vref;
Vrip = discrete_frequency_real_PI.sim_out.Ripple;
Tref = discrete_frequency_real_PI.variable_reference_step_period;
[discrete_R, discrete_RV] = map_voltage_frequency(Vref, Vrip, Tref);

Vref = classic_ripple.sim_out.Vref;
Vrip = classic_ripple.sim_out.Ripple;
Tref = classic_ripple.variable_reference_step_period;
[classic_R,classic_RV] = map_voltage_frequency(Vref, Vrip, Tref);

figure
hold all
plot(continuous_1_RV, continuous_1_R*100./continuous_1_RV, 'LineWidth',1.5,'color','blue', 'DisplayName', 'RNS Controller')
plot(continuous_RV, continuous_R*100./continuous_RV, 'LineWidth',1.5,'color',continuous_color, 'DisplayName', 'QNS Controller')
plot(discrete_RV, discrete_R*100./discrete_RV, 'LineWidth',1.5,'color',discrete_color, 'DisplayName', 'RS Controller')
plot(classic_RV, classic_R*100./classic_RV, 'LineWidth',1.5,'color',classic_color, 'DisplayName', 'PI Controller')
% plot(continuous_RV, smooth(continuous_R))
% plot(discrete_RV, smooth(discrete_R))
hold off
le = legend;
le.FontSize = 15;
%set(legend, 'Position',[0.345238102475802 0.866269843540494 0.333928564190865 0.0869047596341088]);
ylabel('Voltage ripple (V)')
xlabel('Output voltage (V)')
%set(gcf, 'Position',  [100, 100, 1600, 780])
ax = gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;
set(gcf,'renderer','Painters')
saveas(gcf, strcat(image_folder, 'ripple_simulation'), image_format);


%% Reference Step - Output Voltage - Continuous Controller 1

continuous_Vout = continuous_1_reference_step.sim_out.Vout;
classic_Vout = classic_reference_step.sim_out.Vout;
experiment_Vout = exp_step_con1_100.Vout;

experiment_Vout.Data = experiment_Vout.Data - 0.34;

continuous_Vout.Time = continuous_Vout.Time - reference_start_time;
classic_Vout.Time = classic_Vout.Time - reference_start_time;

figure
hold on
plot(continuous_Vout, 'LineWidth',1.5,'color','blue', 'DisplayName', 'RNS Simulation')
plot(experiment_Vout, 'LineWidth',1.5,'color','red', 'DisplayName', 'RNS Experiment')
plot(classic_Vout, '--', 'LineWidth',1.5,'color','black', 'DisplayName', 'PI SImulation')
hold off
le = legend;
le.FontSize = 15;
set(le, 'Position',[0.602380983343436 0.392460319730971 0.199999996753676 0.0869047596341087]);
ylabel('Output voltage (V)')
xlabel('Time (s)')
ylim([0 110])
xlim([-0.005 0.09])
grid
title('')
ax = gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;
set(gcf,'renderer','Painters')
saveas(gcf, strcat(image_folder, 'step_continuous_1'), image_format);


%% Reference Step - Output Voltage - Continuous Controller 2

continuous_Vout = continuous_2_reference_step.sim_out.Vout;
classic_Vout = classic_reference_step.sim_out.Vout;
experiment_Vout = exp_step_con2_100.Vout;

experiment_Vout.Data = experiment_Vout.Data + 0.13;

continuous_Vout.Time = continuous_Vout.Time - reference_start_time;
classic_Vout.Time = classic_Vout.Time - reference_start_time;

figure
hold on
plot(continuous_Vout, 'LineWidth',1.5,'color','blue', 'DisplayName', 'QNS Simulation')
plot(experiment_Vout, 'LineWidth',1.5,'color','red', 'DisplayName', 'QNS Experiment')
plot(classic_Vout, '--', 'LineWidth',1.5,'color','black', 'DisplayName', 'PI Simulation')
hold off
le = legend;
le.FontSize = 15;
set(le, 'Position',[0.602380983343436 0.392460319730971 0.199999996753676 0.0869047596341087]);
ylabel('Output voltage (V)')
xlabel('Time (s)')
ylim([0 110])
xlim([-0.005 0.45])
grid
ax = gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;
set(gcf,'renderer','Painters')
saveas(gcf, strcat(image_folder, 'step_continuous_2'), image_format);


%% Reference Step - Output Voltage - Discrete Controller

discrete_Vout = discrete_reference_step.sim_out.Vout;
experiment_Vout = exp_step_dis_100.Vout;
classic_Vout = classic_reference_step.sim_out.Vout;

experiment_Vout.Data = experiment_Vout.Data + 0.13;

discrete_Vout.Time = discrete_Vout.Time - reference_start_time;
classic_Vout.Time = classic_Vout.Time - reference_start_time;

figure
hold on
plot(discrete_Vout, 'LineWidth',1.5,'color','blue', 'DisplayName', 'RS Simulation')
plot(experiment_Vout, 'LineWidth',1.5,'color','red', 'DisplayName', 'RS Experiment')
plot(classic_Vout, '--', 'LineWidth',1.5,'color','black', 'DisplayName', 'PI Simulation')
hold off
le = legend;
le.FontSize = 15;
set(le, 'Position',[0.602380983343436 0.392460319730971 0.199999996753676 0.0869047596341087]);
ylabel('Output voltage (V)')
xlabel('Time (s)')
ylim([0 110])
xlim([-0.005 0.09])
grid
ax = gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;
set(gcf,'renderer','Painters')
saveas(gcf, strcat(image_folder, 'step_discrete'), image_format);


%% Reference Step - Output Voltage - Classic Controller

classic_Vout = classic_reference_step.sim_out.Vout;
experiment_Vout = exp_step_pwm_100.Vout;

experiment_Vout.Data = experiment_Vout.Data + 0.13;

classic_Vout.Time = classic_Vout.Time - reference_start_time;

figure
hold on
plot(classic_Vout, 'LineWidth',1.5,'color','blue', 'DisplayName', 'Simulation')
plot(experiment_Vout, 'LineWidth',1.5,'color','red', 'DisplayName', 'Experiment')
hold off
le = legend;
le.FontSize = 15;
set(le, 'Position',[0.602380983343436 0.392460319730971 0.199999996753676 0.0869047596341087]);
ylabel('Output voltage (V)')
xlabel('Time (s)')
ylim([0 110])
xlim([0 0.16])
title('')
grid
ax = gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;
set(gcf,'renderer','Painters')
saveas(gcf, strcat(image_folder, 'step_classic'), image_format);


%% Reference Step - Inductor Current - Continuous Controller 1

continuous_IL = continuous_1_reference_step.sim_out.IL;
experiment_IL = exp_step_con1_100.IL;

experiment_IL.Data = experiment_IL.Data + 0.35;

continuous_IL.Time = continuous_IL.Time - reference_start_time;
experiment_IL.Time = experiment_IL.Time;

figure
hold on
plot(continuous_IL, 'LineWidth',1.5,'color','blue', 'DisplayName', 'Simulation')
plot(experiment_IL, 'LineWidth',1.5,'color','red', 'DisplayName', 'Experiment')
hold off
le = legend;
le.FontSize = 15;
ylabel('Coil Current (A)')
xlabel('Time (s)')
xlim([-0.005 0.09])
ylim([0 25])
title('')
grid
ax = gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;
set(gcf,'renderer','Painters')
saveas(gcf, strcat(image_folder, 'step_current_continuous_1'), image_format);


%% Reference Step - Inductor Current - Continuous Controller 2

continuous_IL = continuous_2_reference_step.sim_out.IL;
experiment_IL = exp_step_con2_100.IL;

experiment_IL.Data = experiment_IL.Data + 0.35;

experiment_IL.Time = experiment_IL.Time - 9.5237e-04;

continuous_IL.Time = continuous_IL.Time - reference_start_time;

figure
hold on
plot(continuous_IL, 'LineWidth',1.5,'color','blue', 'DisplayName', 'Simulation')
plot(experiment_IL, 'LineWidth',1.5,'color','red', 'DisplayName', 'Experiment')
hold off
le = legend;
le.FontSize = 15;
ylabel('Coil Current (A)')
xlabel('Time (s)')
xlim([-0.005 0.45])
ylim([0 10])
grid
ax = gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;
set(gcf,'renderer','Painters')
saveas(gcf, strcat(image_folder, 'step_current_continuous_2'), image_format);


%% Reference Step - Inductor Current - Discrete Controller

discrete_IL = discrete_reference_step.sim_out.IL;
experiment_IL = exp_step_dis_100.IL;

experiment_IL.Data = experiment_IL.Data + 0.35;

experiment_IL.Time = experiment_IL.Time - 6.8430e-04;
discrete_IL.Time = discrete_IL.Time - reference_start_time;

figure
hold on
plot(discrete_IL, 'LineWidth',1.5,'color','blue', 'DisplayName', 'Simulation')
plot(experiment_IL, 'LineWidth',1.5,'color','red', 'DisplayName', 'Experiment')
hold off
le = legend;
le.FontSize = 15;
ylabel('Coil Current (A)')
xlabel('Time (s)')
ylim([0 25])
xlim([-0.005 0.09])
grid
ax = gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;
set(gcf,'renderer','Painters')
saveas(gcf, strcat(image_folder, 'step_current_discrete'), image_format);


%% Reference Step - Inductor Current - Classic Controller

classic_IL = classic_reference_step.sim_out.IL;
experiment_IL = exp_step_pwm_100.IL;

experiment_IL.Data = experiment_IL.Data + 0.43;

classic_IL.Time = classic_IL.Time - reference_start_time;
experiment_IL.Time = experiment_IL.Time + 5e-5;

figure
hold on
plot(classic_IL, 'LineWidth',1.5,'color','blue', 'DisplayName', 'Simulation')
plot(experiment_IL, 'LineWidth',1.5,'color','red', 'DisplayName', 'Experiment')
hold off
le = legend;
le.FontSize = 15;
ylabel('Coil Current (A)')
xlabel('Time (s)')
ylim([0 25])
xlim([0 0.18])
title('')
grid
ax = gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;
set(gcf,'renderer','Painters')
saveas(gcf, strcat(image_folder, 'step_current_classic'), image_format);


%% Load Step - Continuous 1

classic_Vout = classic_load_step.sim_out.Vout;
continuous_Vout = continuous_1_load_step.sim_out.Vout;
experiment_Vout = exp_load_con1_C_100.Vout;

classic_Vout.Time = classic_Vout.Time - classic_load_step.disturbance_Ro_time;
continuous_Vout.Time = continuous_Vout.Time - continuous_1_load_step.disturbance_Ro_time;

experiment_Vout.Data = experiment_Vout.Data - 1.5;
experiment_Vout.Time = experiment_Vout.Time - 0.02;

Vref = classic_load_step.test_voltages;

figure
hold on
plot(experiment_Vout, 'DisplayName', 'RNS Experiment', 'color', 'red')
plot(classic_Vout, 'DisplayName', 'PI Simulation', 'color', 'black')
plot(continuous_Vout, 'DisplayName', 'RNS Simulation', 'color', 'blue')
le = legend;
le.FontSize = 15;
hold off
ylim([Vref*0.96 Vref*1.04])
xlim([-0.042 0.158])
ylabel('Output voltage (V)')
xlabel('Time (s)')
grid
ax = gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;
set(gcf,'renderer','Painters')
saveas(gcf, strcat(image_folder, 'load_step_correction_continous_1'), image_format);


%% Load Step - Continuous 2

classic_Vout = classic_load_step.sim_out.Vout;
continuous_Vout = continuous_2_load_step.sim_out.Vout;
experiment_Vout = exp_load_con2_C_100.Vout;

classic_Vout.Time = classic_Vout.Time - classic_load_step.disturbance_Ro_time;
continuous_Vout.Time = continuous_Vout.Time - continuous_2_load_step.disturbance_Ro_time;

experiment_Vout.Data = experiment_Vout.Data - 1.5;
experiment_Vout.Time = experiment_Vout.Time - 0.022;

Vref = classic_load_step.test_voltages;

figure
hold on
plot(experiment_Vout, 'DisplayName', 'QNS Experiment', 'color', 'red')
plot(classic_Vout, 'DisplayName', 'PI Simulation', 'color', 'black')
plot(continuous_Vout, 'DisplayName', 'QNS Simulation', 'color', 'blue')
le = legend;
le.FontSize = 15;
hold off
ylim([Vref*0.96 Vref*1.04])
xlim([-0.042 0.158])
ylabel('Output voltage (V)')
xlabel('Time (s)')
grid
ax = gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;
set(gcf,'renderer','Painters')
saveas(gcf, strcat(image_folder, 'load_step_correction_continous_2'), image_format);


%% Load Step - Discrete

classic_Vout = classic_load_step.sim_out.Vout;
discrete_Vout = discrete_load_step.sim_out.Vout;
experiment_Vout = exp_load_dis_C_100.Vout;

classic_Vout.Time = classic_Vout.Time - classic_load_step.disturbance_Ro_time;
discrete_Vout.Time = discrete_Vout.Time - discrete_load_step.disturbance_Ro_time;

experiment_Vout.Data = experiment_Vout.Data - 1.5;
experiment_Vout.Time = experiment_Vout.Time - 0.02;

Vref = classic_load_step.test_voltages;

figure
hold on
plot(experiment_Vout, 'DisplayName', 'RS Experiment', 'color', 'red')
plot(classic_Vout, 'DisplayName', 'PI Simulation', 'color', 'black')
plot(discrete_Vout, 'DisplayName', 'RS Simulation', 'color', 'blue')
le = legend;
le.FontSize = 15;
hold off
ylim([Vref*0.96 Vref*1.04])
xlim([-0.04 0.158])
ylabel('Output voltage (V)')
xlabel('Time (s)')
grid
ax = gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;
set(gcf,'renderer','Painters')
saveas(gcf, strcat(image_folder, 'load_step_correction_discrete'), image_format);

%%

step_QNS = get_step_info(continuous_2_step_analysis.sim_out);
step_RNS = get_step_info(continuous_1_step_analysis.sim_out);
step_RS = get_step_info(discrete_step_analysis.sim_out);
step_PWM = get_step_info(classic_step_analysis.sim_out);


figure
hold all
plot(step_QNS.Vref, [step_QNS.Vinfo.SettlingTime], 'LineWidth',1.5, 'DisplayName', 'QNS Controller')
plot(step_RNS.Vref, [step_RNS.Vinfo.SettlingTime], 'LineWidth',1.5, 'DisplayName', 'RNS Controller')
plot(step_RS.Vref, [step_RS.Vinfo.SettlingTime], 'LineWidth',1.5, 'DisplayName', 'RS Controller')
plot(step_PWM.Vref, [step_PWM.Vinfo.SettlingTime], 'LineWidth',1.5, 'DisplayName', 'PWM Controller')
hold off
le = legend;
le.FontSize = 15;
le.Position = [0.172023817080827 0.651984121510231 0.34464284958584 0.260714278334663];
ylabel('Settling Time (s)')
xlabel('Voltage reference (V)')
grid
ax = gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;
set(gcf,'renderer','Painters')
saveas(gcf, strcat(image_folder, 'step_settling_time'), image_format);


figure
hold all
plot(step_QNS.Vref, [step_QNS.Iinfo.Peak], 'LineWidth',1.5, 'DisplayName', 'QNS Controller')
plot(step_RNS.Vref, [step_RNS.Iinfo.Peak], 'LineWidth',1.5, 'DisplayName', 'RNS Controller')
plot(step_RS.Vref, [step_RS.Iinfo.Peak], 'LineWidth',1.5, 'DisplayName', 'RS Controller')
plot(step_PWM.Vref, [step_PWM.Iinfo.Peak], 'LineWidth',1.5, 'DisplayName', 'PWM Controller')
hold off
le = legend;
le.FontSize = 15;
le.Position = [0.172023817080827 0.651984121510231 0.34464284958584 0.260714278334663];
ylabel('Current Peak (A)')
xlabel('Voltage reference (V)')
grid
ax = gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;
set(gcf,'renderer','Painters')
saveas(gcf, strcat(image_folder, 'step_peak_current'), image_format);

%%

step_QNS = get_step_info(continuous_2_step_analysis.sim_out);
step_RNS = get_step_info(continuous_1_step_analysis.sim_out);
step_RS = get_step_info(discrete_step_analysis.sim_out);
step_PWM = get_step_info(classic_step_analysis.sim_out);


figure
hold all
plot(step_QNS.Vref, [step_QNS.Vinfo.SettlingTime], 'LineWidth',1.5, 'DisplayName', 'QNS Controller')
plot(step_RNS.Vref, [step_RNS.Vinfo.SettlingTime], 'LineWidth',1.5, 'DisplayName', 'RNS Controller')
plot(step_RS.Vref, [step_RS.Vinfo.SettlingTime], 'LineWidth',1.5, 'DisplayName', 'RS Controller')
plot(step_PWM.Vref, [step_PWM.Vinfo.SettlingTime], 'LineWidth',1.5, 'DisplayName', 'PWM Controller')
hold off
le = legend;
le.FontSize = 12;
le.Position = [0.172023817080827 0.651984121510231 0.34464284958584 0.260714278334663];
ylabel('Settling Time (s)')
xlabel('Voltage reference (V)')
grid
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
set(gcf,'renderer','Painters')
saveas(gcf, strcat(image_folder, 'step_settling_time'), image_format);


figure
hold all
plot(step_QNS.Vref, [step_QNS.Iinfo.Peak], 'LineWidth',1.5, 'DisplayName', 'QNS Controller')
plot(step_RNS.Vref, [step_RNS.Iinfo.Peak], 'LineWidth',1.5, 'DisplayName', 'RNS Controller')
plot(step_RS.Vref, [step_RS.Iinfo.Peak], 'LineWidth',1.5, 'DisplayName', 'RS Controller')
plot(step_PWM.Vref, [step_PWM.Iinfo.Peak], 'LineWidth',1.5, 'DisplayName', 'PWM Controller')
hold off
le = legend;
le.FontSize = 12;
le.Position = [0.172023817080827 0.651984121510231 0.34464284958584 0.260714278334663];
ylabel('Current Peak (A)')
xlabel('Voltage reference (V)')
grid
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
set(gcf,'renderer','Painters')
saveas(gcf, strcat(image_folder, 'step_peak_current'), image_format);


%%
data = {};

Vref = continuous_1_frequency_ideal_40.sim_out.Vref;
Verr = continuous_1_frequency_ideal_40.sim_out.Verr;
Tref = continuous_1_frequency_ideal_40.variable_reference_step_period;
[continuous_1_40_Er, continuous_1_40_Vr] = map_voltage_frequency(Vref, Verr, Tref);

Vref = continuous_2_frequency_ideal_40.sim_out.Vref;
Verr = continuous_2_frequency_ideal_40.sim_out.Verr;
Tref = continuous_2_frequency_ideal_40.variable_reference_step_period;
[continuous_2_40_Er, continuous_2_40_Vr] = map_voltage_frequency(Vref, Verr, Tref);

Vref = discrete_frequency_ideal_40.sim_out.Vref;
Verr = discrete_frequency_ideal_40.sim_out.Verr;
Tref = discrete_frequency_ideal_40.variable_reference_step_period;
[discrete_40_Er, discrete_40_Vr] = map_voltage_frequency(Vref, Verr, Tref);

Vref = continuous_1_frequency_ideal_200.sim_out.Vref;
Verr = continuous_1_frequency_ideal_200.sim_out.Verr;
Tref = continuous_1_frequency_ideal_200.variable_reference_step_period;
[continuous_1_200_Er, continuous_1_200_Vr] = map_voltage_frequency(Vref, Verr, Tref);

Vref = continuous_2_frequency_ideal_200.sim_out.Vref;
Verr = continuous_2_frequency_ideal_200.sim_out.Verr;
Tref = continuous_2_frequency_ideal_200.variable_reference_step_period;
[continuous_2_200_Er, continuous_2_200_Vr] = map_voltage_frequency(Vref, Verr, Tref);

Vref = discrete_frequency_ideal_200.sim_out.Vref;
Verr = discrete_frequency_ideal_200.sim_out.Verr;
Tref = discrete_frequency_ideal_200.variable_reference_step_period;
[discrete_200_Er, discrete_200_Vr] = map_voltage_frequency(Vref, Verr, Tref);

Vref = continuous_1_frequency_ideal_1M.sim_out.Vref;
Verr = continuous_1_frequency_ideal_1M.sim_out.Verr;
Tref = continuous_1_frequency_ideal_1M.variable_reference_step_period;
[continuous_1_1M_Er, continuous_1_1M_Vr] = map_voltage_frequency(Vref, Verr, Tref);

Vref = continuous_2_frequency_ideal_1M.sim_out.Vref;
Verr = continuous_2_frequency_ideal_1M.sim_out.Verr;
Tref = continuous_2_frequency_ideal_1M.variable_reference_step_period;
[continuous_2_1M_Er, continuous_2_1M_Vr] = map_voltage_frequency(Vref, Verr, Tref);

Vref = discrete_frequency_ideal_1M.sim_out.Vref;
Verr = discrete_frequency_ideal_1M.sim_out.Verr;
Tref = discrete_frequency_ideal_1M.variable_reference_step_period;
[discrete_1M_Er, discrete_1M_Vr] = map_voltage_frequency(Vref, Verr, Tref);

data.c1_40 = mean(abs(100*continuous_1_40_Er./continuous_1_40_Vr));
data.c2_40 = mean(abs(100*continuous_2_40_Er./continuous_2_40_Vr));
data.d1_40 = mean(abs(100*discrete_40_Er./discrete_40_Vr));

data.c1_200 = mean(abs(100*continuous_1_200_Er./continuous_1_200_Vr));
data.c2_200 = mean(abs(100*continuous_2_200_Er./continuous_2_200_Vr));
data.d1_200 = mean(abs(100*discrete_200_Er./discrete_200_Vr));

data.c1_1M = mean(abs(100*continuous_1_1M_Er./continuous_1_1M_Vr));
data.c2_1M = mean(abs(100*continuous_2_1M_Er./continuous_2_1M_Vr));
data.d1_1M = mean(abs(100*discrete_1M_Er./discrete_1M_Vr));

printErrorFrequencyTable(data)
