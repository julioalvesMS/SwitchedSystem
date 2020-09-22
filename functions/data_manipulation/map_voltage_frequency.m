function [F,FV] = map_voltage_frequency(Vref,Freq,Tref)

    if (~exist('Tref','var'))
        Tref = 0.5;
    end
    
    dt = Vref.TimeInfo.Increment;
    dt = Freq.Time(end) - Freq.Time(end-1);
    F = [];
    FV = [];
    
    for t=0.99:Tref:max(Freq.Time)
        data = getsampleusingtime(Freq, t-Tref*0.2, t+dt);
        Rdata = getsampleusingtime(Vref, t-Tref*0.2, t+dt);

        F(end+1) = mean(data.Data);
        FV(end+1) = mean(Rdata.Data);
    end
end

