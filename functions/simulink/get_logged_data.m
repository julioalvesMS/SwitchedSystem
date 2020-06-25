function var = get_logged_data(logsout, var_name, compression_rate)

    if ~exist('compression_rate', 'var')
        compression_rate = 1;
    end
    
    if ~isempty(logsout.find(var_name))
        var = downsample_timeseries(logsout.get(var_name).Values, compression_rate);
    else
        var = timeseries(var_name);
    end
    
end

