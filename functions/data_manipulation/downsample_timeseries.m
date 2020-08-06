function out = downsample_timeseries(data, compression_factor)

    out = resample(data, downsample(data.Time, compression_factor));
    if isnan(out.TimeInfo.Increment)
        interval = out.Time(end) - out.Time(end-1);
        out = setuniformtime(out,'Interval',interval);
    end
end

