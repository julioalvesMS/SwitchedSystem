function plot_voltage_lambda_equilibrium(equilibrium_points, lambda, name, folder)
%PLOT_VOLTAGE_CURRENT_EQUILIBRIUM Plot Voltage x Current
%   Makes a Voltage X Current plot from the simulation output.
%   The image will be saved automatically if the folder is argument 
%   is given

    configuration = figure_configuration;
    
    configuration.save = nargin == 4;
    
    if configuration.save
        configuration.folder_path = folder;
    end
    
    configuration.title = strcat(name, ' - Saída x Duty Cycle');

    data.x = lambda(:,1);
    data.y = equilibrium_points(:,2);

    
    configuration.ylabel = 'v_o [V]';
    configuration.xlabel = '\lambda_1';
    
    plot_figure(data, configuration);
end

