function plot_voltage_time(sim_out, name, folder)
%PLOT_VOLTAGE_CURRENT Plot Voltage x Current
%   Makes a Voltage X Current plot from the simulation output.
%   The image will be saved automatically if the folder is argument 
%   is given
    
    configuration = figure_configuration;
    
    configuration.save = nargin == 3;
    
    if configuration.save
        configuration.folder_path = folder;
    end
    
    configuration.title = strcat(name, ' - Tensão x Tempo');

    for i=length(sim_out):-1:1
        data(2*i).x = sim_out(i).Vout.Time*1e3;
        data(2*i).y = sim_out(i).Vout.Data;
        
        data(2*i-1).x = sim_out(i).xe.Time*1e3;
        data(2*i-1).y = sim_out(i).xe.Data(:,2);
    end
    
    configuration.ylabel = 'v_o [V]';
    configuration.xlabel = 't [ms]';
    
    plot_figure(data, configuration);
end

