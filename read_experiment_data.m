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
N = length(data_files);

for i = 3:N
    file = data_files(i);
    load(file.name);
end

folder = "CURVAS_2020_01_14";

folder_path = fullfile(root_data_folder, folder);

data_files = dir(folder_path);
N = length(data_files);

for i = 3:N
    file = data_files(i);
    load(file.name);
end


%% Process useful data

ratio = 1;

% Chaveamento - Contínuo 1

C1_BuckBoost = {};
C1_BuckBoost.t = C1BuckBoost_C1_noRef_Step_30V_Vo00000(round(9*end/20):ratio:end,1);
C1_BuckBoost.Vo = C1BuckBoost_C1_noRef_Step_30V_Vo00000(round(9*end/20):ratio:end,2);
C1_BuckBoost.Vref = C3BuckBoost_C1_noRef_Step_30V_Vref00000(round(9*end/20):ratio:end,2);
C1_BuckBoost.IL = C2BuckBoost_C1_noRef_Step_30V_IL00000(round(9*end/20):ratio:end,2);

C1_Boost = {};
C1_Boost.t = C1Boost_C1_noRef_Step_90V_Vo00000(round(9*end/20):ratio:end,1);
C1_Boost.Vo = C1Boost_C1_noRef_Step_90V_Vo00000(round(9*end/20):ratio:end,2);
C1_Boost.Vref = C3Boost_C1_noRef_Step_90V_Vref00000(round(9*end/20):ratio:end,2);
C1_Boost.IL = C2Boost_C1_noRef_Step_90V_IL00000(round(9*end/20):ratio:end,2);

C1_Buck = {};
C1_Buck.t = C1Buck_C1_noRef_Step_30V_Vo00000(round(9*end/20):ratio:end,1);
C1_Buck.Vo = C1Buck_C1_noRef_Step_30V_Vo00000(round(9*end/20):ratio:end,2);
C1_Buck.Vref = C3Buck_C1_noRef_Step_30V_Vref00000(round(9*end/20):ratio:end,2);
C1_Buck.IL = C2Buck_C1_noRef_Step_30V_IL00000(round(9*end/20):ratio:end,2);


% Chaveamento - Contínuo 2


C2_BuckBoost = {};
C2_BuckBoost.t = C1BuckBoost_C2_noRef_Step_30V_Vo00000(round(9*end/20):ratio:end,1);
C2_BuckBoost.Vo = C1BuckBoost_C2_noRef_Step_30V_Vo00000(round(9*end/20):ratio:end,2);
C2_BuckBoost.Vref = C3BuckBoost_C2_noRef_Step_30V_Vref00000(round(9*end/20):ratio:end,2);
C2_BuckBoost.IL = C2BuckBoost_C2_noRef_Step_30V_IL00000(round(9*end/20):ratio:end,2);

C2_Boost = {};
C2_Boost.t = C1Boost_C2_noRef_Step_90V_Vo00000(round(9*end/20):ratio:end,1);
C2_Boost.Vo = C1Boost_C2_noRef_Step_90V_Vo00000(round(9*end/20):ratio:end,2);
C2_Boost.Vref = C3Boost_C2_noRef_Step_90V_Vref00000(round(9*end/20):ratio:end,2);
C2_Boost.IL = C2Boost_C2_noRef_Step_90V_IL00000(round(9*end/20):ratio:end,2);

C2_Buck = {};
C2_Buck.t = C1Buck_C2_noRef_Step_30V_Vo00000(round(9*end/20):ratio:end,1);
C2_Buck.Vo = C1Buck_C2_noRef_Step_30V_Vo00000(round(9*end/20):ratio:end,2);
C2_Buck.Vref = C3Buck_C2_noRef_Step_30V_Vref00000(round(9*end/20):ratio:end,2);
C2_Buck.IL = C2Buck_C2_noRef_Step_30V_IL00000(round(9*end/20):ratio:end,2);


% Chaveamento - Contínuo 2 com ref


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


% Chaveamento - PWM


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


% Buck Boost 3


C1_BuckBoost3 = {};
C1_BuckBoost3.t = C1BuckBoost3_C1_noRef_Step_30V_Vo00000(round(9*end/20):ratio:end,1);
C1_BuckBoost3.Vo = C1BuckBoost3_C1_noRef_Step_30V_Vo00000(round(9*end/20):ratio:end,2);
C1_BuckBoost3.Vref = C3BuckBoost3_C1_noRef_Step_30V_Vref00000(round(9*end/20):ratio:end,2);
C1_BuckBoost3.IL = C2BuckBoost3_C1_noRef_Step_30V_IL00000(round(9*end/20):ratio:end,2);

C1R_BuckBoost3 = {};
C1R_BuckBoost3.t = C1BuckBoost3_C1_Vo_Step_30V00000(round(9*end/20):ratio:end,1);
C1R_BuckBoost3.Vo = C1BuckBoost3_C1_Vo_Step_30V00000(round(9*end/20):ratio:end,2);
C1R_BuckBoost3.Vref = C3BuckBoost3_C1_Step_30V_Vref00000(round(9*end/20):ratio:end,2);
C1R_BuckBoost3.IL = C2BuckBoost3_C1_Step_30V_IL00000(round(9*end/20):ratio:end,2);

C2_BuckBoost3 = {};
C2_BuckBoost3.t = C1BuckBoost3_C2_noRef_Step_30V_Vo00001(round(9*end/20):ratio:end,1);
C2_BuckBoost3.Vo = C1BuckBoost3_C2_noRef_Step_30V_Vo00001(round(9*end/20):ratio:end,2);
C2_BuckBoost3.IL = C4BuckBoost3_C2_noRef_Step_30V_IL00001(round(9*end/20):ratio:end,2);

C2R_BuckBoost3 = {};
C2R_BuckBoost3.t = C1BuckBoost3_C2_Step_30V_Vo00000(round(9*end/20):ratio:end,1);
C2R_BuckBoost3.Vo = C1BuckBoost3_C2_Step_30V_Vo00000(round(9*end/20):ratio:end,2);
C2R_BuckBoost3.Vref = C3BuckBoost3_C2_Step_30V_Vref00000(round(9*end/20):ratio:end,2);
C2R_BuckBoost3.IL = C2BuckBoost3_C2_Step_30V_IL00000(round(9*end/20):ratio:end,2);

%% Remove noise

PWM_Buck.Vof = smooth(PWM_Buck.t,PWM_Buck.Vo,10,'rloess');
PWM_Boost.Vof = smooth(PWM_Boost.t,PWM_Boost.Vo,10,'rloess');
PWM_BuckBoost.Vof = smooth(PWM_BuckBoost.t,PWM_BuckBoost.Vo,10,'rloess');
C2R_Buck.Vof = smooth(C2R_Buck.t,C2R_Buck.Vo,10,'rloess');
C2R_Boost.Vof = smooth(C2R_Boost.t,C2R_Boost.Vo,10,'rloess');
C2R_BuckBoost.Vof = smooth(C2R_BuckBoost.t,C2R_BuckBoost.Vo,10,'rloess');
C2_Buck.Vof = smooth(C2_Buck.t,C2_Buck.Vo,10,'rloess');
C2_Boost.Vof = smooth(C2_Boost.t,C2_Boost.Vo,10,'rloess');
C2_BuckBoost.Vof = smooth(C2_BuckBoost.t,C2_BuckBoost.Vo,10,'rloess');
C1_Buck.Vof = smooth(C1_Buck.t,C1_Buck.Vo,10,'rloess');
C1_Boost.Vof = smooth(C1_Boost.t,C1_Boost.Vo,10,'rloess');
C1_BuckBoost.Vof = smooth(C1_BuckBoost.t,C1_BuckBoost.Vo,10,'rloess');


C1_BuckBoost3.Vof = smooth(C1_BuckBoost3.t,C1_BuckBoost3.Vo,10,'rloess');
C2_BuckBoost3.Vof = smooth(C2_BuckBoost3.t,C2_BuckBoost3.Vo,10,'rloess');
C1R_BuckBoost3.Vof = smooth(C1R_BuckBoost3.t,C1R_BuckBoost3.Vo,10,'rloess');
C2R_BuckBoost3.Vof = smooth(C2R_BuckBoost3.t,C2R_BuckBoost3.Vo,10,'rloess');

%%

% close all
% figure
% hold on
% plot(C2R_BuckBoost.Vo)
% plot(smooth(C2R_BuckBoost.t,C2R_BuckBoost.Vo,10,'lowess'))
% hold off
% title('lowess')
% 
% figure
% hold on
% plot(C2R_BuckBoost.Vo)
% plot(smooth(C2R_BuckBoost.t,C2R_BuckBoost.Vo,10,'sgolay'))
% hold off
% title('sgolay')
% 
% figure
% hold on
% plot(PWM_Boost.Vo)
% plot(smooth(PWM_Boost.t,PWM_Boost.Vo,10,'rloess'))
% hold off
% title('rloess')


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


% C2

% Voltage
plot_experiment_voltage(C2_Buck, 30)
set(gcf,'renderer','Painters')
saveas(gcf, strcat(root_image_folder, 'Experiment_C2_Buck_Vout'), 'epsc');
        
plot_experiment_voltage(C2_Boost, 90)
set(gcf,'renderer','Painters')
saveas(gcf, strcat(root_image_folder, 'Experiment_C2_Boost_Vout'), 'epsc');

plot_experiment_voltage(C2_BuckBoost, 30)
set(gcf,'renderer','Painters')
saveas(gcf, strcat(root_image_folder, 'Experiment_C2_BuckBoost_Vout'), 'epsc');

% Current
plot_experiment_current(C2_Buck)
set(gcf,'renderer','Painters')
saveas(gcf, strcat(root_image_folder, 'Experiment_C2_Buck_IL'), 'epsc');
        
plot_experiment_current(C2_Boost)
set(gcf,'renderer','Painters')
saveas(gcf, strcat(root_image_folder, 'Experiment_C2_Boost_IL'), 'epsc');

plot_experiment_current(C2_BuckBoost)
set(gcf,'renderer','Painters')
saveas(gcf, strcat(root_image_folder, 'Experiment_C2_BuckBoost_IL'), 'epsc');


% C1

% Voltage
plot_experiment_voltage(C1_Buck, 30)
set(gcf,'renderer','Painters')
saveas(gcf, strcat(root_image_folder, 'Experiment_C1_Buck_Vout'), 'epsc');
        
plot_experiment_voltage(C1_Boost, 90)
set(gcf,'renderer','Painters')
saveas(gcf, strcat(root_image_folder, 'Experiment_C1_Boost_Vout'), 'epsc');

plot_experiment_voltage(C1_BuckBoost, 30)
set(gcf,'renderer','Painters')
saveas(gcf, strcat(root_image_folder, 'Experiment_C1_BuckBoost_Vout'), 'epsc');

% Current
plot_experiment_current(C1_Buck)
set(gcf,'renderer','Painters')
saveas(gcf, strcat(root_image_folder, 'Experiment_C1_Buck_IL'), 'epsc');
        
plot_experiment_current(C1_Boost)
set(gcf,'renderer','Painters')
saveas(gcf, strcat(root_image_folder, 'Experiment_C1_Boost_IL'), 'epsc');

plot_experiment_current(C1_BuckBoost)
set(gcf,'renderer','Painters')
saveas(gcf, strcat(root_image_folder, 'Experiment_C1_BuckBoost_IL'), 'epsc');


% Buck Boost 3

% Voltage
plot_experiment_voltage(C1_BuckBoost3, 30)
set(gcf,'renderer','Painters')
saveas(gcf, strcat(root_image_folder, 'Experiment_C1_BuckBoost3_Vout'), 'epsc');

plot_experiment_voltage(C2_BuckBoost3, 30)
set(gcf,'renderer','Painters')
saveas(gcf, strcat(root_image_folder, 'Experiment_C2_BuckBoost3_Vout'), 'epsc');
        
plot_experiment_voltage(C1R_BuckBoost3, 30)
set(gcf,'renderer','Painters')
saveas(gcf, strcat(root_image_folder, 'Experiment_C1R_BuckBoost3_Vout'), 'epsc');

plot_experiment_voltage(C2R_BuckBoost3, 30)
set(gcf,'renderer','Painters')
saveas(gcf, strcat(root_image_folder, 'Experiment_C2R_BuckBoost3_Vout'), 'epsc');

% Current
plot_experiment_current(C1_BuckBoost3)
set(gcf,'renderer','Painters')
saveas(gcf, strcat(root_image_folder, 'Experiment_C1_BuckBoost3_IL'), 'epsc');

plot_experiment_current(C2_BuckBoost3)
set(gcf,'renderer','Painters')
saveas(gcf, strcat(root_image_folder, 'Experiment_C2_BuckBoost3_IL'), 'epsc');
        
plot_experiment_current(C1R_BuckBoost3)
set(gcf,'renderer','Painters')
saveas(gcf, strcat(root_image_folder, 'Experiment_C1R_BuckBoost3_IL'), 'epsc');

plot_experiment_current(C2R_BuckBoost3)
set(gcf,'renderer','Painters')
saveas(gcf, strcat(root_image_folder, 'Experiment_C2R_BuckBoost3_IL'), 'epsc');


