function plot_figure(data, configuration)
%PLOT_FIGURE Summary of this function goes here
%   Detailed explanation goes here
    

    figure;
    hold all;
    
    for i=1:length(data)
        plot(data(i).x1, data(i).y1);
        if isfield(data, 'y2')
            plot(data(i).x2, data(i).y2, 'black--');
        end
    end
    
    ylabel(configuration.ylabel);
    xlabel(configuration.xlabel);
    
    if ~isempty(configuration.title)
        title(configuration.title);
    end
    if ~isempty(configuration.legend)
        legend(configuration.legend);
    end
    
    hold off;
    
    if configuration.save
        
        if isempty(configuration.file_name)
            configuration.file_name = remove_special_characters(configuration.title);
        end
        
        set(gcf,'renderer','Painters')
        saveas(gcf, strcat(configuration.folder_path, configuration.file_name), 'epsc');
    end
end

