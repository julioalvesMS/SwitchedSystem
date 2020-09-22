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
circuit_buck_boost_3_stage = buck_boost_non_inverting(R, Ro, Co, L);


%% Run PWM

converter_template_folder = fullfile('CCS', 'templates', 'Converter');

config_buck.circuit = circuit_buck;
config_buck.template = fullfile(converter_template_folder, 'buck.cpp');
config_boost.circuit = circuit_boost;
config_boost.template = fullfile(converter_template_folder, 'boost.cpp');
config_buck_boost.circuit = circuit_buck_boost;
config_buck_boost.template = fullfile(converter_template_folder, 'buck_boost.cpp');
config_buck_boost_3_stage.circuit = circuit_buck_boost_3_stage;
config_buck_boost_3_stage.template = fullfile(converter_template_folder, 'buck_boost_3.cpp');

% config_converters = {config_buck, config_boost, config_buck_boost, config_buck_boost_3_stage};
config_converters = {config_buck_boost};

%%


bar = waitbar(0, 'Inicializing', 'name', 'Code Generation');

Ns = length(config_converters);
try
    for conv_id=1:Ns
        config = config_converters{conv_id};
        circuit = config.circuit;

        waitbar((conv_id-1)/Ns, bar, sprintf('Converter %d of %d\n%s',conv_id, Ns, circuit.name));

        run load_circuit_sys

        range = circuit.operation_range_voltage_min:5:circuit.operation_range_voltage_max;
        lambdas = generate_lambda_voltage(sys, range);
        dlambdas = generate_lambda_voltage(dsys, range);

        % Calculate the P matrix, as the equilibrium point. Calculation will be
        % in accordance with the chosen control theorem

        [Pc1] = calc_sys_theorem_1_range(sys, lambdas);
        [Pc2] = calc_sys_theorem_2(sys);
        [Pd1, gamma] = calc_sys_discrete_theorem_1_test_1(sys, dsys, dlambdas);

        
        Vref = circuit.limit_cycle_voltage;
        xe = calculate_equilibrium_point(circuit, Vs, Vref);
%         xe = [0; Vref];
        Gamma = circuit.limit_cycle_gamma;
        [cand, kappa] = find_cycles(dsys, xe, Gamma);

        cycle_cost = find_limit_cycle(dsys, kappa, cand, 1);
        cycle_H2 = find_limit_cycle(dsys, kappa, cand, 2);
        cycle_Hinf = find_limit_cycle(dsys, kappa, cand, 3);
        
        
        A = [-R/L  -1/L
              1/Co  -1/(Ro*Co)
        ];
        B = [Vs/L; 0];
        C = [1.8 1];
        D = 0;
        sys_ss = ss(A, B, C, D);
        [K,M] = project_state_feedback_h2(sys_ss, pwm_period);


        getP.fnc = gen_fnc_getP(circuit.class_name, Pc1, Pc2, Pd1);
%         getH.fnc = gen_fnc_getH(circuit.class_name, hd1);
%         getD.fnc = gen_fnc_getD(circuit.class_name, dd1);
        DefineSystem.fnc = gen_fnc_DefineSystem(sys);
        DefineDiscreteSystem.fnc = gen_fnc_DefineDiscreteSystem(dsys);
        getClassicVoltageController.fnc = gen_fnc_getClassicVoltageController(circuit, pwm_period);
        getClassicVoltageCurrentController.fnc = gen_fnc_getClassicVoltageCurrentController(circuit, pwm_period);
        getStateFeedbackH2Controller.fnc = gen_fnc_getStateFeedbackH2Controller(circuit, K, C, M);
        getReferenceController.fnc = gen_fnc_getReferenceController(circuit, Tref);
        getCurrentCorrectionController.fnc = gen_fnc_getCurrentCorrectionController(circuit, Tref);
        DefineLimitCycleCost.fnc = gen_fnc_DefineLimitCycleCost(dsys, cycle_cost);
        DefineLimitCycleH2.fnc = gen_fnc_DefineLimitCycleH2(dsys, cycle_H2);
        DefineLimitCycleHinf.fnc = gen_fnc_DefineLimitCycleHinf(dsys, cycle_Hinf);

        getP.tag = 'getP';
        DefineSystem.tag = 'DefineSystem';
        DefineDiscreteSystem.tag = 'DefineDiscreteSystem';
        getClassicVoltageController.tag = 'getClassicVoltageController';
        getClassicVoltageCurrentController.tag = 'getClassicVoltageCurrentController';
        getStateFeedbackH2Controller.tag = 'getStateFeedbackH2Controller';
        getReferenceController.tag = 'getReferenceController';
        getCurrentCorrectionController.tag = 'getCurrentCorrectionController';
        DefineLimitCycleCost.tag = 'DefineLimitCycleCost';
        DefineLimitCycleH2.tag = 'DefineLimitCycleH2';
        DefineLimitCycleHinf.tag = 'DefineLimitCycleHinf';

        file_fnc = {
            getP
            getClassicVoltageController
            getClassicVoltageCurrentController
            getStateFeedbackH2Controller
            getReferenceController
            getCurrentCorrectionController
            DefineSystem
            DefineDiscreteSystem
            DefineLimitCycleCost
            DefineLimitCycleH2
            DefineLimitCycleHinf
            };

        template = fileread(config.template);

        out = template;
        for i=1:length(file_fnc)
            tag = strcat('{{', file_fnc{i}.tag, '}}');
            out = regexprep(out, tag, file_fnc{i}.fnc);
        end

        template_file = dir(config.template);
        file_name = strcat(root_output_folder, template_file.name);
        f = fopen(file_name, 'w');
        fwrite(f, out);
        fclose(f);
    end
catch exception
    close(bar);
    rethrow(exception);
end
close(bar);
