function data = get_step_info(sim_out)
    data = {};
    for i=length(sim_out):-1:1
        Vout = sim_out(i).Vout;
        IL = sim_out(i).IL;
        Vref = sim_out(i).Vref;
        data.Vinfo(i) = stepinfo(Vout.Data, Vout.Time);
        data.Iinfo(i) = stepinfo(IL.Data, IL.Time);
        data.Vref(i) = Vref.Data(end);
    end
end

