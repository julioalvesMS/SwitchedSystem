%% Initial Setup
clear; clc; close all;

% folders to create
root_output_folder = 'CCS/output_files';
cache_folder = 'tmp/cache';

[~,~]=mkdir(root_output_folder);
root_output_folder = strcat(root_output_folder, '/');
[~,~]=mkdir(cache_folder);
cache_folder = strcat(cache_folder, '/');

addpath(genpath('functions'))
addpath(genpath('simulations'))
addpath(genpath('system_models'))
addpath(genpath('scripts'))
addpath(genpath('CCS'))

Simulink.fileGenControl('set', 'CacheFolder', cache_folder);

%% System Specifications

run system_specifications

%% Default Parameters

% Desired DC-DC converter to use
% Options can be found in the system directory:
%   buck
%   boost
%   buck_boost
circuit_buck = buck(R, Ro, Co, L);
circuit_boost = boost(R, Ro, Co, L);
circuit_buck_boost = buck_boost(R, Ro, Co, L);


%% Run PWM


config_converters = {circuit_buck, circuit_boost, circuit_buck_boost};


%%
for conv_id=1:length(config_converters)


    circuit = config_converters{conv_id};
    run load_circuit_sys
    
    lambdas = generate_lambda_voltage(sys, circuit.single_voltage);

    % Calculate the P matrix, as the equilibrium point. Calculation will be
    % in accordance with the chosen control theorem
    
	[Pc1] = calc_sys_theorem_1(sys, lambdas);
    
	[Pc2] = calc_sys_theorem_2(sys, lambdas);
    
	[Pd1, hd1, dd1, ~, dsys] = calc_sys_discrete_theorem_1(sys, dsys, lambdas);
    
    
    fnc_getP = gen_fnc_getP(circuit.class_name, Pc1, Pc2, Pd1);
    fnc_getH = gen_fnc_getH(circuit.class_name, hd1);
    fnc_getD = gen_fnc_getD(circuit.class_name, dd1);
    fnc_DefineDiscreteSystem = gen_fnc_DefineDiscreteSystem(dsys);
    
    file_fnc = {
        fnc_getP
        fnc_getH
        fnc_getD
        fnc_DefineDiscreteSystem
        };
    
    str=sprintf('%s\n\n',file_fnc{:});
    
    file_name = strcat(root_output_folder, circuit.class_name, ".cpp");
    f = fopen(file_name, 'w');
    fwrite(f, str);
    fclose(f);
end




