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

folder = "CURVAS_2020_02_13";

folder_path = fullfile(root_data_folder, folder, '*.dat');

data_files = dir(folder_path);
N = length(data_files);

for i = 1:N
    file = data_files(i);
    load(file.name);
end


%% Process useful data

ratio = 1;

% Transitorio
Buck_H2 = {};
Buck_H2.t = C1Buck_H2_T00000(1:ratio:end,1);
Buck_H2.TST = C1Buck_H2_T00000(1:ratio:end,2);
Buck_H2.SW = C2Buck_H2_T00000(1:ratio:end,2);
Buck_H2.Vo = C3Buck_H2_T00000(1:ratio:end,2);
Buck_H2.IL = C4Buck_H2_T00000(1:ratio:end,2);

Buck_PWM = {};
Buck_PWM.t = C1Buck_PWM_T00000(1:ratio:end,1);
Buck_PWM.TST = C1Buck_PWM_T00000(1:ratio:end,2);
Buck_PWM.SW = C2Buck_PWM_T00000(1:ratio:end,2);
Buck_PWM.Vo = C3Buck_PWM_T00000(1:ratio:end,2);
Buck_PWM.IL = C4Buck_PWM_T00000(1:ratio:end,2);
%% Remove noise


Buck_H2.Vof = smooth(Buck_H2.t,Buck_H2.Vo,10,'rloess');
Buck_PWM.Vof = smooth(Buck_PWM.t,Buck_PWM.Vo,10,'rloess');


%% Plot


figure
hold on
plot(Buck_PWM.t, Buck_PWM.Vo)
plot(Buck_PWM.t, Buck_PWM.Vof)


%%


t = Buck_PWM.t;
l = ones(size(t));
figure
hold on
plot(Buck_PWM.t, Buck_PWM.Vof)
plot(Buck_CL_H2_T.t, Buck_CL_H2_T.Vof)
plot(Buck_CA_T.t, Buck_CA_T.Vof)
plot(t, l*30*1.02, '--b')
plot(t, l*30*0.98, '--b')
legend({'PI', 'CL H2', 'CA'})


figure
hold on
plot(Buck_PWM.t, Buck_PWM.IL)
plot(Buck_CL_H2_T.t, Buck_CL_H2_T.IL)
plot(Buck_CA_T.t, Buck_CA_T.IL)
legend({'PI', 'CL H2', 'CA'})

%%


t = Buck_PI_R.t;
l = ones(size(t));
figure
hold on
plot(Buck_PI_R.t, Buck_PI_R.Vof)
plot(Buck_CL_H2_R.t, Buck_CL_H2_R.Vof)
plot(Buck_CA_R.t, Buck_CA_R.Vof)
plot(t, l*30*1.02, '--b')
plot(t, l*30*0.98, '--b')
legend({'PI', 'CL H2', 'CA'})


figure
hold on
plot(Buck_PI_R.t, Buck_PI_R.IL)
plot(Buck_CL_H2_R.t, Buck_CL_H2_R.IL)
plot(Buck_CA_R.t, Buck_CA_R.IL)
legend({'PI', 'CL H2', 'CA'})