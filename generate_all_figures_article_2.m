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
data_folder = "CURVAS_2020_07_27";

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

file = 'all_simulations_article_2.mat';

load(file)

%% Load experiment data

folder_path = fullfile(root_data_folder, data_folder, '*.dat');

data_files = dir(folder_path);
N = length(data_files);
 
for i = 1:N
    file = data_files(i);
    load(file.name);
end

ratio = 1;

% Degrau - Clássico - 100V
data = {};
data_C1 = C1STEP_PWM_100V00000;
data_C4 = C4STEP_PWM_100V00000;
data.t = data_C1(1:ratio:end,1);
data.Vout = timeseries(data_C1(1:ratio:end,2), data.t);
data.Vref = timeseries(100*ones(size(data.Vout)), data.t);
data.IL = timeseries(data_C4(1:ratio:end,2), data.t);
exp_step_pwm_100 = data;

% Degrau - Contínuo 1 - 100V
data = {};
data_C1 = C1STEP_CON1_100V00000;
data_C4 = C4STEP_CON1_100V00000;
data.t = data_C1(1:ratio:end,1);
data.Vout = timeseries(data_C1(1:ratio:end,2), data.t);
data.Vref = timeseries(100*ones(size(data.Vout)), data.t);
data.IL = timeseries(data_C4(1:ratio:end,2), data.t);
exp_step_con1_100 = data;

% Degrau - Contínuo 2 - 100V
data = {};
data_C1 = C1STEP_CON2_100V00002;
data_C4 = C4STEP_CON2_100V00002;
data.t = data_C1(1:ratio:end,1);
data.Vout = timeseries(data_C1(1:ratio:end,2), data.t);
data.Vref = timeseries(100*ones(size(data.Vout)), data.t);
data.IL = timeseries(data_C4(1:ratio:end,2), data.t);
exp_step_con2_100 = data;

% Degrau - Discreto - 100V
data = {};
data_C1 = C1STEP_DIS_100V00001;
data_C4 = C4STEP_DIS_100V00001;
data.t = data_C1(1:ratio:end,1);
data.Vout = timeseries(data_C1(1:ratio:end,2), data.t);
data.Vref = timeseries(100*ones(size(data.Vout)), data.t);
data.IL = timeseries(data_C4(1:ratio:end,2), data.t);
exp_step_dis_100 = data;




% Degrau Carga - Clássico - 100V
data = {};
data_C1 = C1LOAD_PWM_100V00000;
data_C4 = C4LOAD_PWM_100V00000;
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
data_C1 = C1LOAD_CON1_C_100V00000;
data_C4 = C4LOAD_CON1_C_100V00000;
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
data_C1 = C1LOAD_CON2_C_100V00000;
data_C4 = C4LOAD_CON2_C_100V00000;
data.t = data_C1(1:ratio:end,1);
data.Vout = timeseries(data_C1(1:ratio:end,2), data.t);
data.Vref = timeseries(100*ones(size(data.Vout)), data.t);
data.IL = timeseries(data_C4(1:ratio:end,2), data.t);
exp_load_con2_C_100 = data;

% Degrau Carga - Discreto - sem Correção de corente - 100V
data = {};
data_C1 = C1LOAD_DIS_noC_100V00101;
data_C4 = C4LOAD_DIS_noC_100V00101;
data.t = data_C1(1:ratio:end,1);
data.Vout = timeseries(data_C1(1:ratio:end,2), data.t);
data.Vref = timeseries(100*ones(size(data.Vout)), data.t);
data.IL = timeseries(data_C4(1:ratio:end,2), data.t);
exp_load_dis_noC_100 = data;

% Degrau Carga - Discreto - com Correção de corente - 100V
data = {};
data_C1 = C1LOAD_DIS_C_100V00101;
data_C4 = C4LOAD_DIS_C_100V00101;
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
dis = csvread('C:\Users\Julio-LEPO\Codes\data\exp_frequency_disc_MH_C.csv');
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
plot(equilibrium(:,1), equilibrium(:,2), 'black','LineWidth',2)
ylabel('v_o [V]');
xlabel('i_L [A]');
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
plot(continuous_2_FV, continuous_2_F/1e3, 'DisplayName', 'Robust Continuous', 'LineWidth',1)
plot(continuous_1_FV, continuous_1_F/1e3, 'DisplayName', 'Restrict Continuous', 'LineWidth',1)
plot(discrete_FV, discrete_F/1e3, 'DisplayName', 'Discrete', 'LineWidth',1)
hold off
ylabel('Switching frequency [kHz]')
xlabel('Output voltage [V]')
legend
% legend1 = legend;
% set(legend1, 'Position',[0.514880959618659 0.176984130390106 0.333928564190864 0.165476185934884]);
%set(gcf, 'Position',  [100, 100, 1600, 780])
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
plot(continuous_2_FV, continuous_2_F/1e3, 'DisplayName', 'Robust Continuous', 'LineWidth',1)
plot(continuous_1_FV, continuous_1_F/1e3, 'DisplayName', 'Restrict Continuous', 'LineWidth',1)
plot(discrete_FV, discrete_F/1e3, 'DisplayName', 'Discrete', 'LineWidth',1)
hold off
ylabel('Switching frequency [kHz]')
xlabel('Output voltage [V]')
legend
% legend1 = legend;
% set(legend1, 'Position',[0.514880959618659 0.176984130390106 0.333928564190864 0.165476185934884]);
%set(gcf, 'Position',  [100, 100, 1600, 780])
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
plot(continuous_2_FV, continuous_2_F/1e3, 'DisplayName', 'Robust Continuous', 'LineWidth',1)
plot(continuous_1_FV, continuous_1_F/1e3, 'DisplayName', 'Restrict Continuous', 'LineWidth',1)
plot(discrete_FV, discrete_F/1e3, 'DisplayName', 'Discrete', 'LineWidth',1)
hold off
ylabel('Switching frequency [kHz]')
xlabel('Output voltage [V]')
legend
% legend1 = legend;
% set(legend1, 'Position',[0.514880959618659 0.176984130390106 0.333928564190864 0.165476185934884]);
%set(gcf, 'Position',  [100, 100, 1600, 780])
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
plot(continuous_2_Vr, continuous_2_Er, 'DisplayName', 'Robust Continuous', 'LineWidth',1)
plot(continuous_1_Vr, continuous_1_Er, 'DisplayName', 'Restrict Continuous', 'LineWidth',1)
plot(discrete_Vr, discrete_Er, 'DisplayName', 'Discrete', 'LineWidth',1)
hold off
grid
ylabel('Voltage Error [V]')
xlabel('Output voltage [V]')
legend
% legend1 = legend;
% set(legend1, 'Position',[0.514880959618659 0.176984130390106 0.333928564190864 0.165476185934884]);
%set(gcf, 'Position',  [100, 100, 1600, 780])
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
plot(continuous_2_Vr, continuous_2_Er, 'DisplayName', 'Robust Continuous', 'LineWidth',1)
plot(continuous_1_Vr, continuous_1_Er, 'DisplayName', 'Restrict Continuous', 'LineWidth',1)
plot(discrete_Vr, discrete_Er, 'DisplayName', 'Discrete', 'LineWidth',1)
hold off
grid
ylabel('Voltage Error [V]')
xlabel('Output voltage [V]')
legend
% legend1 = legend;
% set(legend1, 'Position',[0.514880959618659 0.176984130390106 0.333928564190864 0.165476185934884]);
%set(gcf, 'Position',  [100, 100, 1600, 780])
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
plot(continuous_2_Vr, continuous_2_Er*100./continuous_2_Vr, 'DisplayName', 'Robust Continuous', 'LineWidth',1)
plot(continuous_1_Vr, continuous_1_Er*100./continuous_2_Vr, 'DisplayName', 'Restrict Continuous', 'LineWidth',1)
plot(discrete_Vr, discrete_Er*100./continuous_2_Vr, 'DisplayName', 'Discrete', 'LineWidth',1)
hold off
grid
ylabel('Voltage Error [V]')
xlabel('Output voltage [V]')
legend
% legend1 = legend;
% set(legend1, 'Position',[0.514880959618659 0.176984130390106 0.333928564190864 0.165476185934884]);
%set(gcf, 'Position',  [100, 100, 1600, 780])
set(gcf,'renderer','Painters')
saveas(gcf, strcat(image_folder, 'voltage_error_40kHz'), image_format);


%% Frequency - Real 40kHz - PI - Continuous Controller 1

Vref = continuous_1_frequency_ideal_PI_40.sim_out.Vref;
Freq = continuous_1_frequency_ideal_PI_40.sim_out.F;
Tref = continuous_1_frequency_ideal_PI_40.variable_reference_step_period;
[switching_ideal_F, switching_ideal_FV] = map_voltage_frequency(Vref, Freq, Tref);

Vref = continuous_1_frequency_real_PI.sim_out.Vref;
Freq = continuous_1_frequency_real_PI.sim_out.F;
Tref = continuous_1_frequency_real_PI.variable_reference_step_period;
[switching_real_F, switching_real_FV] = map_voltage_frequency(Vref, Freq, Tref);

switching_classic_FV = switching_real_FV;
switching_classic_F = 20e3*ones(size((switching_classic_FV)));

exp_FV = exp_frequency_MH_C.Vref;
exp_F = exp_frequency_noMH_noC.Fcon1;

figure
hold all
plot(switching_ideal_FV, switching_ideal_F/1e3, '--', 'DisplayName', 'Switched: Ideal Switches', 'LineWidth',1,'color','blue')
plot(switching_real_FV, switching_real_F/1e3, '-x', 'DisplayName', 'Switched: Real Switches', 'LineWidth',1,'color','blue')
plot(switching_classic_FV, switching_classic_F/1e3, '--', 'DisplayName', 'PWM', 'LineWidth',1,'color','black')
plot(exp_FV, exp_F, '-o', 'DisplayName', 'Experimental Result', 'color', 'red')
hold off
legend1 = legend;
set(legend1, 'Position',[0.514880959618659 0.176984130390106 0.333928564190864 0.165476185934884]);
ylabel('Switching frequency [kHz]')
xlabel('Output voltage [V]')
%set(gcf, 'Position',  [100, 100, 1600, 780])
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
switching_classic_F = 20e3*ones(size((switching_classic_FV)));

exp_FV = exp_frequency_MH_C.Vref;
exp_F = exp_frequency_MH_C.Fcon2;

figure
hold all
plot(switching_ideal_FV, switching_ideal_F/1e3, '--', 'DisplayName', 'Switched: Ideal Switches', 'LineWidth',1,'color','blue')
plot(switching_real_FV, switching_real_F/1e3, '-x', 'DisplayName', 'Switched: Real Switches', 'LineWidth',1,'color','blue')
plot(switching_classic_FV, switching_classic_F/1e3, '--', 'DisplayName', 'PWM', 'LineWidth',1,'color','black')
plot(exp_FV, exp_F, '-o', 'DisplayName', 'Experimental Result', 'color', 'red')
hold off
legend1 = legend;
set(legend1, 'Position',[0.514880959618659 0.176984130390106 0.333928564190864 0.165476185934884]);
ylabel('Switching frequency [kHz]')
xlabel('Output voltage [V]')
%set(gcf, 'Position',  [100, 100, 1600, 780])
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
switching_classic_F = 20e3*ones(size((switching_classic_FV)));

exp_FV = exp_frequency_MH_C.Vref;
exp_F = exp_frequency_MH_C.Fdis;

figure
hold all
plot(switching_ideal_FV, switching_ideal_F/1e3, '--', 'DisplayName', 'Switched: Ideal Switches', 'LineWidth',1,'color','blue')
plot(switching_real_FV, switching_real_F/1e3, '-x', 'DisplayName', 'Switched: Real Switches', 'LineWidth',1,'color','blue')
plot(switching_classic_FV, switching_classic_F/1e3, '--', 'DisplayName', 'PWM', 'LineWidth',1,'color','black')
plot(exp_FV, exp_F, '-o', 'DisplayName', 'Experimental Result', 'color', 'red')
hold off
legend1 = legend;
set(legend1, 'Position',[0.514880959618659 0.176984130390106 0.333928564190864 0.165476185934884]);
ylabel('Switching frequency [kHz]')
xlabel('Output voltage [V]')
%set(gcf, 'Position',  [100, 100, 1600, 780])
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
plot(continuous_1_RV, continuous_1_R*100./continuous_1_RV, 'LineWidth',1.5,'color','blue', 'DisplayName', 'Continuous 1')
plot(continuous_RV, continuous_R*100./continuous_RV, 'LineWidth',1.5,'color',continuous_color, 'DisplayName', 'Continuous 2')
plot(discrete_RV, discrete_R*100./discrete_RV, 'LineWidth',1.5,'color',discrete_color, 'DisplayName', 'Discrete')
plot(classic_RV, classic_R*100./classic_RV, 'LineWidth',1.5,'color',classic_color, 'DisplayName', 'Classic')
% plot(continuous_RV, smooth(continuous_R))
% plot(discrete_RV, smooth(discrete_R))
hold off
legend
%set(legend, 'Position',[0.345238102475802 0.866269843540494 0.333928564190865 0.0869047596341088]);
ylabel('Voltage Ripple [V]')
xlabel('Output voltage [V]')
%set(gcf, 'Position',  [100, 100, 1600, 780])
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
plot(continuous_Vout, 'LineWidth',2,'color','blue', 'DisplayName', 'Simulation')
plot(experiment_Vout, 'LineWidth',2,'color','red', 'DisplayName', 'Experiment')
plot(classic_Vout, '--', 'LineWidth',2,'color','black', 'DisplayName', 'Classic')
hold off
legend1 = legend;
set(legend1, 'Position',[0.602380983343436 0.392460319730971 0.199999996753676 0.0869047596341087]);
ylabel('V_{out} [V]')
xlabel('t [s]')
ylim([0 110])
xlim([-0.005 0.09])
grid
title('')
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
plot(continuous_Vout, 'LineWidth',2,'color','blue', 'DisplayName', 'Simulation')
plot(experiment_Vout, 'LineWidth',2,'color','red', 'DisplayName', 'Experiment')
plot(classic_Vout, '--', 'LineWidth',2,'color','black', 'DisplayName', 'Classic')
hold off
legend1 = legend;
set(legend1, 'Position',[0.602380983343436 0.392460319730971 0.199999996753676 0.0869047596341087]);
ylabel('V_{out} [V]')
xlabel('t [s]')
ylim([0 110])
xlim([-0.005 0.45])
grid
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
plot(discrete_Vout, 'LineWidth',2,'color','blue', 'DisplayName', 'Simulation')
plot(experiment_Vout, 'LineWidth',2,'color','red', 'DisplayName', 'Experiment')
plot(classic_Vout, '--', 'LineWidth',2,'color','black', 'DisplayName', 'PWM')
hold off
legend1 = legend;
set(legend1, 'Position',[0.602380983343436 0.392460319730971 0.199999996753676 0.0869047596341087]);
ylabel('V_{out} [V]')
xlabel('t [s]')
ylim([0 110])
xlim([-0.005 0.09])
grid
set(gcf,'renderer','Painters')
saveas(gcf, strcat(image_folder, 'step_discrete'), image_format);


%% Reference Step - Output Voltage - Classic Controller

classic_Vout = classic_reference_step.sim_out.Vout;
experiment_Vout = exp_step_pwm_100.Vout;

experiment_Vout.Data = experiment_Vout.Data + 0.13;

classic_Vout.Time = classic_Vout.Time - reference_start_time;

figure
hold on
plot(classic_Vout, 'LineWidth',2,'color','blue', 'DisplayName', 'Simulation')
plot(experiment_Vout, 'LineWidth',2,'color','red', 'DisplayName', 'Experiment')
hold off
legend1 = legend;
set(legend1, 'Position',[0.602380983343436 0.392460319730971 0.199999996753676 0.0869047596341087]);
ylabel('V_{out} [V]')
xlabel('t [s]')
ylim([0 110])
xlim([0 0.16])
title('')
grid
set(gcf,'renderer','Painters')
saveas(gcf, strcat(image_folder, 'step_classic'), image_format);


%% Reference Step - Inductor Current - Continuous Controller 1

continuous_IL = continuous_1_reference_step.sim_out.IL;
experiment_IL = exp_step_con1_100.IL;

experiment_IL.Data = experiment_IL.Data + 0.3;

continuous_IL.Time = continuous_IL.Time - reference_start_time;
experiment_IL.Time = experiment_IL.Time - 0.000358;

figure
hold on
plot(continuous_IL, 'LineWidth',2,'color','blue', 'DisplayName', 'Simulation')
plot(experiment_IL, 'LineWidth',2,'color','red', 'DisplayName', 'Experiment')
hold off
legend;
ylabel('Coil Current [A]')
xlabel('t [s]')
xlim([-0.005 0.09])
ylim([0 25])
title('')
grid
set(gcf,'renderer','Painters')
saveas(gcf, strcat(image_folder, 'step_current_continuous_1'), image_format);


%% Reference Step - Inductor Current - Continuous Controller 2

continuous_IL = continuous_2_reference_step.sim_out.IL;
experiment_IL = exp_step_con2_100.IL;

experiment_IL.Data = experiment_IL.Data + 0.3;

continuous_IL.Time = continuous_IL.Time - reference_start_time;

figure
hold on
plot(continuous_IL, 'LineWidth',2,'color','blue', 'DisplayName', 'Simulation')
plot(experiment_IL, 'LineWidth',2,'color','red', 'DisplayName', 'Experiment')
hold off
legend1 = legend;
ylabel('Coil Current [A]')
xlabel('t [s]')
xlim([-0.005 0.18])
ylim([0 10])
grid
set(gcf,'renderer','Painters')
saveas(gcf, strcat(image_folder, 'step_current_continuous_2'), image_format);


%% Reference Step - Inductor Current - Discrete Controller

discrete_IL = discrete_reference_step.sim_out.IL;
experiment_IL = exp_step_dis_100.IL;

experiment_IL.Data = experiment_IL.Data + 0.37;

discrete_IL.Time = discrete_IL.Time - reference_start_time;

figure
hold on
plot(discrete_IL, 'LineWidth',2,'color','blue', 'DisplayName', 'Simulation')
plot(experiment_IL, 'LineWidth',2,'color','red', 'DisplayName', 'Experiment')
hold off
legend
ylabel('Coil Current [A]')
xlabel('t [s]')
ylim([0 20])
xlim([-0.005 0.09])
grid
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
plot(classic_IL, 'LineWidth',2,'color','blue', 'DisplayName', 'Simulation')
plot(experiment_IL, 'LineWidth',2,'color','red', 'DisplayName', 'Experiment')
hold off
legend
ylabel('Coil Current [A]')
xlabel('t [s]')
ylim([0 20])
xlim([0 0.18])
title('')
grid
set(gcf,'renderer','Painters')
saveas(gcf, strcat(image_folder, 'step_current_classic'), image_format);


%% Load Step - Continuous 1

classic_Vout = classic_load_step.sim_out.Vout;
continuous_Vout = continuous_1_load_step.sim_out.Vout;
experiment_Vout = exp_load_con1_C_100.Vout;

classic_Vout.Time = classic_Vout.Time - classic_load_step.disturbance_Ro_time;
continuous_Vout.Time = continuous_Vout.Time - continuous_1_load_step.disturbance_Ro_time;

experiment_Vout.Data = experiment_Vout.Data - 0.9;
experiment_Vout.Time = experiment_Vout.Time - 0.022;

Vref = classic_load_step.test_voltages;

figure
hold on
plot(experiment_Vout, 'DisplayName', 'Experimental Result', 'color', 'red')
plot(classic_Vout, 'DisplayName', 'Classic Controller', 'color', 'black')
plot(continuous_Vout, 'DisplayName', 'Retrict Continuous Switched Controller', 'color', 'blue')
legend
hold off
ylim([Vref*0.92 Vref*1.06])
xlim([-0.042 0.158])
ylabel('V_{out} [V]')
xlabel('t [s]')
grid
set(gcf,'renderer','Painters')
saveas(gcf, strcat(image_folder, 'load_step_correction_continous_1'), image_format);


%% Load Step - Continuous 2

classic_Vout = classic_load_step.sim_out.Vout;
continuous_Vout = continuous_2_load_step.sim_out.Vout;
experiment_Vout = exp_load_con2_C_100.Vout;

classic_Vout.Time = classic_Vout.Time - classic_load_step.disturbance_Ro_time;
continuous_Vout.Time = continuous_Vout.Time - continuous_2_load_step.disturbance_Ro_time;

experiment_Vout.Data = experiment_Vout.Data - 0.9;
experiment_Vout.Time = experiment_Vout.Time - 0.022;

Vref = classic_load_step.test_voltages;

figure
hold on
plot(experiment_Vout, 'DisplayName', 'Experimental Result', 'color', 'red')
plot(classic_Vout, 'DisplayName', 'Classic Controller', 'color', 'black')
plot(continuous_Vout, 'DisplayName', 'Robust Continuous Switched Controller', 'color', 'blue')
legend
hold off
ylim([Vref*0.92 Vref*1.06])
xlim([-0.042 0.158])
ylabel('V_{out} [V]')
xlabel('t [s]')
grid
set(gcf,'renderer','Painters')
saveas(gcf, strcat(image_folder, 'load_step_correction_continous_2'), image_format);


%% Load Step - Discrete

classic_Vout = classic_load_step.sim_out.Vout;
discrete_Vout = discrete_load_step.sim_out.Vout;
experiment_Vout = exp_load_dis_C_100.Vout;

classic_Vout.Time = classic_Vout.Time - classic_load_step.disturbance_Ro_time;
discrete_Vout.Time = discrete_Vout.Time - discrete_load_step.disturbance_Ro_time;

experiment_Vout.Data = experiment_Vout.Data - 1.9;
experiment_Vout.Time = experiment_Vout.Time - 0.022;

Vref = classic_load_step.test_voltages;

figure
hold on
plot(experiment_Vout, 'DisplayName', 'Experimental Result', 'color', 'red')
plot(classic_Vout, 'DisplayName', 'Classic Controller', 'color', 'black')
plot(discrete_Vout, 'DisplayName', 'Discrete Switched Controller', 'color', 'blue')
legend
hold off
ylim([Vref*0.96 Vref*1.04])
xlim([-0.04 0.158])
ylabel('V_{out} [V]')
xlabel('t [s]')
grid
set(gcf,'renderer','Painters')
saveas(gcf, strcat(image_folder, 'load_step_correction_discrete'), image_format);