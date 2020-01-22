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

folder = "CURVAS_2020_01_22";

folder_path = fullfile(root_data_folder, folder, '*.dat');

data_files = dir(folder_path);
N = length(data_files);

for i = 1:N
    file = data_files(i);
    load(file.name);
end


%% Process useful data

ratio = 1;

% Boost
CO_Buck = {};
CO_Buck.t = C3Buck_LC_CO_40V_TS00000(1:ratio:end,1);
CO_Buck.Vo = C1Buck_LC_CO_40V_Vo00000(1:ratio:end,2);
CO_Buck.SW = C2Buck_LC_CO_40V_SW00000(1:ratio:end,2);
CO_Buck.TS = C3Buck_LC_CO_40V_TS00000(1:ratio:end,2);
CO_Buck.IL = C4Buck_LC_CO_40V_IL00000(1:ratio:end,2)*10;
H2_Buck = {};
H2_Buck.t = C3Buck_LC_H2_40V_TS00000(1:ratio:end,1);
H2_Buck.Vo = C1Buck_LC_H2_40V_Vo00000(1:ratio:end,2);
H2_Buck.SW = C2Buck_LC_H2_40V_SW00000(1:ratio:end,2);
H2_Buck.TS = C3Buck_LC_H2_40V_TS00000(1:ratio:end,2);
H2_Buck.IL = C4Buck_LC_H2_40V_IL00000(1:ratio:end,2)*10;
Hi_Buck = {};
Hi_Buck.t = C3Buck_LC_Hi_40V_TS00000(1:ratio:end,1);
Hi_Buck.Vo = C1Buck_LC_Hi_40V_Vo00000(1:ratio:end,2);
Hi_Buck.SW = C2Buck_LC_Hi_40V_SW00000(1:ratio:end,2);
Hi_Buck.TS = C3Buck_LC_Hi_40V_TS00000(1:ratio:end,2);
Hi_Buck.IL = C4Buck_LC_Hi_40V_IL00000(1:ratio:end,2)*10;

% Boost
CO_Boost = {};
CO_Boost.t = C3Boost_LC_CO_100V_TS00000(1:ratio:end,1);
CO_Boost.Vo = C1Boost_LC_CO_100V_Vo00000(1:ratio:end,2);
CO_Boost.SW = C2Boost_LC_CO_100V_SW00000(1:ratio:end,2);
CO_Boost.TS = C3Boost_LC_CO_100V_TS00000(1:ratio:end,2);
CO_Boost.IL = C4Boost_LC_CO_100V_IL00000(1:ratio:end,2)*10;

% Buck Boost
CO_BuckBoost = {};
CO_BuckBoost.t = C3BuckBoost_LC_CO_80V_TS00000(1:ratio:end,1);
CO_BuckBoost.Vo = C1BuckBoost_LC_CO_80V_Vo00000(1:ratio:end,2);
CO_BuckBoost.SW = C2BuckBoost_LC_CO_80V_SW00000(1:ratio:end,2);
CO_BuckBoost.TS = C3BuckBoost_LC_CO_80V_TS00000(1:ratio:end,2);
CO_BuckBoost.IL = C4BuckBoost_LC_CO_80V_IL00000(1:ratio:end,2)*10;

% Buck Boost 3
CO_BuckBoost3 = {};
CO_BuckBoost3.t = C3BuckBoost3_LC_CO_80V_TS00000(1:ratio:end,1);
CO_BuckBoost3.Vo = C1BuckBoost3_LC_CO_80V_Vo00000(1:ratio:end,2);
CO_BuckBoost3.SW = C2BuckBoost3_LC_CO_80V_SW00000(1:ratio:end,2);
CO_BuckBoost3.TS = C3BuckBoost3_LC_CO_80V_TS00000(1:ratio:end,2);
CO_BuckBoost3.IL = C4BuckBoost3_LC_CO_80V_IL00000(1:ratio:end,2)*10;

%% Remove noise


CO_Buck.Vof = smooth(CO_Buck.t,CO_Buck.Vo,10,'rloess');
H2_Buck.Vof = smooth(H2_Buck.t,H2_Buck.Vo,10,'rloess');
Hi_Buck.Vof = smooth(Hi_Buck.t,Hi_Buck.Vo,10,'rloess');

CO_Boost.Vof = smooth(CO_Boost.t,CO_Boost.Vo,10,'rloess');

CO_BuckBoost.Vof = smooth(CO_BuckBoost.t,CO_BuckBoost.Vo,10,'rloess');

CO_BuckBoost3.Vof = smooth(CO_BuckBoost3.t,CO_BuckBoost3.Vo,10,'rloess');



%% Plot


figure
plot(CO_Buck.t, CO_Buck.Vo)
figure
plot(CO_Buck.t, CO_Buck.Vof)

figure
plot(CO_Buck.IL, CO_Buck.Vo)
figure
plot(CO_Buck.IL, CO_Buck.Vof)



%%


figure
plot(CO_Buck.IL, CO_Buck.Vof)