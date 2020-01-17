function plot_experiment_voltage(data, ref)
    ratio = 1;

    y = data.Vof(1:ratio:end);
    x = 1e3*data.t(1:ratio:end);
    r = ref*ones(size(y));
    figure
    hold on
    plot(x, y)
    plot(x, r, '--black')
    
    ylabel('V_o [V]')
    xlabel('t [ms]')
    
    axis([-100 max(x) max(min(data.Vo),0) ref*1.5])
    hold off
end
