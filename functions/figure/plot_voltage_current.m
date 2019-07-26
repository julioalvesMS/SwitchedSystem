function plot_voltage_current(sim_out, name, folder)
%PLOT_VOLTAGE_CURRENT Plot Voltage x Current
%   Makes a Voltage X Current plot from the simulation output.
%   The image will be saved automatically if the folder is argument 
%   is given

    configuration = figure_configuration;
    
    configuration.save = nargin == 3;
    
    if configuration.save
        configuration.folder_path = folder;
    end
    
    configuration.title = strcat(name, ' - Tens�o x Corrente');

    for i=length(sim_out):-1:1
        data(i).x = sim_out(i).Iout.Data;
        data(i).y = sim_out(i).Vout.Data;
    end
    
    configuration.ylabel = 'v_o [V]';
    configuration.xlabel = 'I [A]';
    
    plot_figure(data, configuration);
end

