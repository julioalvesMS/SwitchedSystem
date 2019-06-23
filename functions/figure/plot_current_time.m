function plot_current_time(sim_out, name, folder)
%PLOT_VOLTAGE_CURRENT Plot Voltage x Current
%   Makes a Voltage X Current plot from the simulation output.
%   The image will be saved automatically if the folder is argument 
%   is given

    configuration = figure_configuration;
    
    configuration.save = nargin == 3;
    
    if configuration.save
        configuration.folder_path = folder;
    end
    
    configuration.title = strcat(name, ' - Corrente x Tempo');

    for i=length(sim_out):-1:1
        data(i).x = sim_out(i).Iout.Time*1e3;
        data(i).y = sim_out(i).Iout.Data;
    end
    
    configuration.ylabel = 'I [A]';
    configuration.xlabel = 't [ms]';
    
    plot_figure(data, configuration);
end

