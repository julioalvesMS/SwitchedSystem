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

folder = "CURVAS_2020_02_04";

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
Buck_PI_T = {};
Buck_PI_T.t = C1Buck_PI_T_30V_Vo00000(1:ratio:end,1);
Buck_PI_T.Vo = C1Buck_PI_T_30V_Vo00000(1:ratio:end,2);
Buck_PI_T.SW = C2Buck_PI_T_30V_SW00000(1:ratio:end,2);
Buck_PI_T.IL = C4Buck_PI_T_30V_IL00000(1:ratio:end,2);
Buck_CL_H2_T = {};
Buck_CL_H2_T.t = C1Buck_CL_H2_T_30V_Vo00000(1:ratio:end,1);
Buck_CL_H2_T.Vo = C1Buck_CL_H2_T_30V_Vo00000(1:ratio:end,2);
Buck_CL_H2_T.SW = C2Buck_CL_H2_T_30V_SW00000(1:ratio:end,2);
Buck_CL_H2_T.IL = C4Buck_CL_H2_T_30V_IL00000(1:ratio:end,2);
Buck_CA_T = {};
Buck_CA_T.t = C1Buck_CA_T_30V_Vo00000(1:ratio:end,1);
Buck_CA_T.Vo = C1Buck_CA_T_30V_Vo00000(1:ratio:end,2);
Buck_CA_T.SW = C2Buck_CA_T_30V_SW00000(1:ratio:end,2);
Buck_CA_T.IL = C4Buck_CA_T_30V_IL00000(1:ratio:end,2);

% Transitorio
Buck_PI_R = {};
Buck_PI_R.t = C1Buck_PI_R_30V_Vo00000(1:ratio:end,1);
Buck_PI_R.Vo = C1Buck_PI_R_30V_Vo00000(1:ratio:end,2);
Buck_PI_R.SW = C2Buck_PI_R_30V_SW00000(1:ratio:end,2);
Buck_PI_R.IL = C4Buck_PI_R_30V_IL00000(1:ratio:end,2);
Buck_CL_H2_R = {};
Buck_CL_H2_R.t = C1Buck_CL_H2_R_30V_Vo00000(1:ratio:end,1);
Buck_CL_H2_R.Vo = C1Buck_CL_H2_R_30V_Vo00000(1:ratio:end,2);
Buck_CL_H2_R.SW = C2Buck_CL_H2_R_30V_SW00000(1:ratio:end,2);
Buck_CL_H2_R.IL = C4Buck_CL_H2_R_30V_IL00000(1:ratio:end,2);
Buck_CA_R = {};
Buck_CA_R.t = C1Buck_CA_R_30V_Vo00000(1:ratio:end,1);
Buck_CA_R.Vo = C1Buck_CA_R_30V_Vo00000(1:ratio:end,2);
Buck_CA_R.SW = C2Buck_CA_R_30V_SW00000(1:ratio:end,2);
Buck_CA_R.IL = C4Buck_CA_R_30V_IL00000(1:ratio:end,2);

%% Remove noise


Buck_PI_T.Vof = smooth(Buck_PI_T.t,Buck_PI_T.Vo,10,'rloess');
Buck_CL_H2_T.Vof = smooth(Buck_CL_H2_T.t,Buck_CL_H2_T.Vo,10,'rloess');
Buck_CA_T.Vof = smooth(Buck_CA_T.t,Buck_CA_T.Vo,10,'rloess');

Buck_PI_R.Vof = smooth(Buck_PI_R.t,Buck_PI_R.Vo,10,'rloess');
Buck_CL_H2_R.Vof = smooth(Buck_CL_H2_R.t,Buck_CL_H2_R.Vo,10,'rloess');
Buck_CA_R.Vof = smooth(Buck_CA_R.t,Buck_CA_R.Vo,10,'rloess');


%% Plot


figure
hold on
plot(Buck_PI_T.t, Buck_PI_T.Vo)
plot(Buck_PI_T.t, Buck_PI_T.Vof)
figure
hold on
plot(Buck_CL_H2_T.t, Buck_CL_H2_T.Vo)
plot(Buck_CL_H2_T.t, Buck_CL_H2_T.Vof)
figure
hold on
plot(Buck_CA_T.t, Buck_CA_T.Vo)
plot(Buck_CA_T.t, Buck_CA_T.Vof)

figure
hold on
plot(Buck_PI_R.t, Buck_PI_R.Vo)
plot(Buck_PI_R.t, Buck_PI_R.Vof)
figure
hold on
plot(Buck_CL_H2_R.t, Buck_CL_H2_R.Vo)
plot(Buck_CL_H2_R.t, Buck_CL_H2_R.Vof)
figure
hold on
plot(Buck_CA_R.t, Buck_CA_R.Vo)
plot(Buck_CA_R.t, Buck_CA_R.Vof)


%%


t = Buck_PI_T.t;
l = ones(size(t));
figure
hold on
plot(Buck_PI_T.t, Buck_PI_T.Vof)
plot(Buck_CL_H2_T.t, Buck_CL_H2_T.Vof)
plot(Buck_CA_T.t, Buck_CA_T.Vof)
plot(t, l*30*1.02, '--b')
plot(t, l*30*0.98, '--b')
legend({'PI', 'CL H2', 'CA'})


figure
hold on
plot(Buck_PI_T.t, Buck_PI_T.IL)
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