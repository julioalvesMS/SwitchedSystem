function plot_experiment_current(data)
    ratio = 1;

    y = data.IL(1:ratio:end);
    x = 1e3*data.t(1:ratio:end);
    figure
    hold on
    plot(x, y)
    
    ylabel('I_L [A]')
    xlabel('t [ms]')
    
    axis([-100 max(x) min(data.IL) max(data.IL)])
    hold off
end