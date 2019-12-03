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


plot_compression_rate = 2e2;


Simulink.fileGenControl('set', 'CacheFolder', cache_folder);

%% Load data

root_data_folder = 'data';

folder = "CURVAS_2019_11_29";

folder_path = fullfile(root_data_folder, folder);

data_files = dir(folder_path);
N = length(data_files)-2;

for i = 3:N
    file = data_files(i);
    load(file.name);
end


%% Process useful data

ratio = 5;

C2R_BuckBoost = {};
C2R_BuckBoost.t = C1BuckBoost_C2_Step_30V_Vo00000(round(9*end/20):ratio:end,1);
C2R_BuckBoost.Vo = C1BuckBoost_C2_Step_30V_Vo00000(round(9*end/20):ratio:end,2);
C2R_BuckBoost.Vref = C3BuckBoost_C2_Step_30V_Vref00000(round(9*end/20):ratio:end,2);
C2R_BuckBoost.IL = C2BuckBoost_C2_Step_30V_IL00000(round(9*end/20):ratio:end,2);

C2R_Boost = {};
C2R_Boost.t = C1Boost_C2_Step_90V_Vo00000(round(9*end/20):ratio:end,1);
C2R_Boost.Vo = C1Boost_C2_Step_90V_Vo00000(round(9*end/20):ratio:end,2);
C2R_Boost.Vref = C3Boost_C2_Step_90V_Vref00000(round(9*end/20):ratio:end,2);
C2R_Boost.IL = C2Boost_C2_Step_90V_IL00000(round(9*end/20):ratio:end,2);

C2R_Buck = {};
C2R_Buck.t = C1Buck_C2_Step_30V_Vo00000(round(9*end/20):ratio:end,1);
C2R_Buck.Vo = C1Buck_C2_Step_30V_Vo00000(round(9*end/20):ratio:end,2);
C2R_Buck.Vref = C3Buck_C2_Step_30V_Vref00000(round(9*end/20):ratio:end,2);
C2R_Buck.IL = C2Buck_C2_Step_30V_IL00000(round(9*end/20):ratio:end,2);


PWM_BuckBoost = {};
PWM_BuckBoost.t = C1BuckBoost_PWM_Step_30V_Vo00000(round(9*end/20):ratio:end,1);
PWM_BuckBoost.Vo = C1BuckBoost_PWM_Step_30V_Vo00000(round(9*end/20):ratio:end,2);
PWM_BuckBoost.Vref = C3BuckBoost_PWM_Step_30V_Vref00000(round(9*end/20):ratio:end,2);
PWM_BuckBoost.IL = C2BuckBoost_PWM_Step_30V_IL00000(round(9*end/20):ratio:end,2);

PWM_Boost = {};
PWM_Boost.t = C1Boost_PWM_Step_90V_Vo00000(round(9*end/20):ratio:end,1);
PWM_Boost.Vo = C1Boost_PWM_Step_90V_Vo00000(round(9*end/20):ratio:end,2);
PWM_Boost.Vref = C3Boost_PWM_Step_90V_Vref00000(round(9*end/20):ratio:end,2);
PWM_Boost.IL = C2Boost_PWM_Step_90V_IL00000(round(9*end/20):ratio:end,2);

PWM_Buck = {};
PWM_Buck.t = C1Buck_PWM_Step_30V_Vo00000(round(9*end/20):ratio:end,1);
PWM_Buck.Vo = C1Buck_PWM_Step_30V_Vo00000(round(9*end/20):ratio:end,2);
PWM_Buck.Vref = C3Buck_PWM_Step_30V_Vref00000(round(9*end/20):ratio:end,2);
PWM_Buck.IL = C2Buck_PWM_Step_30V_IL00000(round(9*end/20):ratio:end,2);

%% Remove noise

%% Plot

% PWM

% Voltage
plot_experiment_voltage(PWM_Buck, 30)
set(gcf,'renderer','Painters')
saveas(gcf, strcat(root_image_folder, 'Experiment_PWM_Buck_Vout'), 'epsc');
        
plot_experiment_voltage(PWM_Boost, 90)
set(gcf,'renderer','Painters')
saveas(gcf, strcat(root_image_folder, 'Experiment_PWM_Boost_Vout'), 'epsc');

plot_experiment_voltage(PWM_BuckBoost, 30)
set(gcf,'renderer','Painters')
saveas(gcf, strcat(root_image_folder, 'Experiment_PWM_BuckBoost_Vout'), 'epsc');

% Current
plot_experiment_current(PWM_Buck)
set(gcf,'renderer','Painters')
saveas(gcf, strcat(root_image_folder, 'Experiment_PWM_Buck_IL'), 'epsc');
        
plot_experiment_current(PWM_Boost)
set(gcf,'renderer','Painters')
saveas(gcf, strcat(root_image_folder, 'Experiment_PWM_Boost_IL'), 'epsc');

plot_experiment_current(PWM_BuckBoost)
set(gcf,'renderer','Painters')
saveas(gcf, strcat(root_image_folder, 'Experiment_PWM_BuckBoost_IL'), 'epsc');


% C2R

% Voltage
plot_experiment_voltage(C2R_Buck, 30)
set(gcf,'renderer','Painters')
saveas(gcf, strcat(root_image_folder, 'Experiment_C2R_Buck_Vout'), 'epsc');
        
plot_experiment_voltage(C2R_Boost, 90)
set(gcf,'renderer','Painters')
saveas(gcf, strcat(root_image_folder, 'Experiment_C2R_Boost_Vout'), 'epsc');

plot_experiment_voltage(C2R_BuckBoost, 30)
set(gcf,'renderer','Painters')
saveas(gcf, strcat(root_image_folder, 'Experiment_C2R_BuckBoost_Vout'), 'epsc');

% Current
plot_experiment_current(C2R_Buck)
set(gcf,'renderer','Painters')
saveas(gcf, strcat(root_image_folder, 'Experiment_C2R_Buck_IL'), 'epsc');
        
plot_experiment_current(C2R_Boost)
set(gcf,'renderer','Painters')
saveas(gcf, strcat(root_image_folder, 'Experiment_C2R_Boost_IL'), 'epsc');

plot_experiment_current(C2R_BuckBoost)
set(gcf,'renderer','Painters')
saveas(gcf, strcat(root_image_folder, 'Experiment_C2R_BuckBoost_IL'), 'epsc');

function plot_experiment_voltage(data, ref)
    figure
    hold on
    plot(data.t*1e3, data.Vo)
    plot(data.t*1e3, ref*ones(size(data.Vo)), '--black')
    
    ylabel('V_o [V]')
    xlabel('t [ms]')
    
    axis([-100 1e3 max(min(data.Vo),0) ref*1.5])
    hold off
end

function plot_experiment_current(data)
    figure
    hold on
    plot(data.t*1e3, data.IL)
    
    ylabel('I_L [A]')
    xlabel('t [ms]')
    
    axis([-100 1e3 max(min(data.IL),0) max(data.IL)])
    hold off
end
