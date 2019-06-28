function plot_disturbance_voltage_time(sim_out, disturbance_time, name, folder)
%PLOT_VOLTAGE_CURRENT Plot Voltage x Current
%   Makes a Voltage X Current plot from the simulation output.
%   The image will be saved automatically if the folder is argument 
%   is given
    
    configuration = figure_configuration;
    
    configuration.save = nargin == 4;
    
    if configuration.save
        configuration.folder_path = folder;
    end
    
    configuration.title = name;
    
    dt = 100e-3; % [s]

    for i=length(sim_out):-1:1
        
        from = find(sim_out(i).Vout.Time > disturbance_time - dt, 1);
        to = find(sim_out(i).Vout.Time > disturbance_time + 3*dt, 1);
        
        if isempty(to)
            to = length(sim_out(i).Vout.Time);
        end
        
        data(i).x = sim_out(i).Vout.Time(1:to-from+1)*1e3;
        data(i).y = sim_out(i).Vout.Data(from:to);
    end
    
    configuration.ylabel = 'v_o [V]';
    configuration.xlabel = 't [ms]';
    
    plot_figure(data, configuration);
end

