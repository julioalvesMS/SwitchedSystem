function plot_voltage_current(sim_out, name, folder)
%PLOT_VOLTAGE_CURRENT Plot Voltage x Current
%   Makes a Voltage X Current plot from the simulation output.
%   The image will be saved automatically if the folder is argument 
%   is given
    
    save_image = nargin == 3;
    
    name = strcat(name, ' - Voltage x Current');

    figure;
    hold all;
    
    for i=1:length(sim_out)
        plot(sim_out(i).x.Data(:,1), sim_out(i).x.Data(:,2));
    end
    
    ylabel('V [V]');
    xlabel('I [A]');
    title(name);
    hold off;
    
    
    if save_image
        saveas(gcf, strcat(folder, name, '.png'));
    end
end

