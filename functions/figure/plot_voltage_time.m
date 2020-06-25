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
    
    configuration.title = strcat(name, ' - Voltage x Time');

    tmp_legends = cell(1,length(sim_out)*2);
    
    for i=length(sim_out):-1:1
        data(i).x1 = sim_out(i).Vout.Time*1e3;
        data(i).y1 = sim_out(i).Vout.Data;
        tmp_legends{2*i-1} = 'Vo';
        if isfield(sim_out(i), 'Vref')
            data(i).x2 = sim_out(i).Vref.Time*1e3;
            data(i).y2 = sim_out(i).Vref.Data;
            tmp_legends{2*i} = 'Vref';
        end
    end
    
%     configuration.legend = tmp_legends;
    
    configuration.ylabel = 'v_o [V]';
    configuration.xlabel = 't [ms]';
    
    configuration.file_name = remove_special_characters(strcat(name, ' - Tensão x Tempo'));

    plot_figure(data, configuration);
end

