function plot_figure(data, configuration)
%PLOT_FIGURE Summary of this function goes here
%   Detailed explanation goes here
    

    figure;
    hold all;
    
    for i=1:length(data)
        plot(data(i).x, data(i).y);
    end
    
    ylabel(configuration.ylabel);
    xlabel(configuration.xlabel);
    title(configuration.title);
    hold off;
    
    
    if configuration.save
        
        if isempty(configuration.file_name)
            configuration.file_name = remove_special_characters(configuration.title);
        end
        
        saveas(gcf, strcat(configuration.folder_path, configuration.file_name, '.png'));
    end
end

