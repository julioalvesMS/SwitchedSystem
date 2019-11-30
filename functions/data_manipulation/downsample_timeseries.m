function out = downsample_timeseries(data, compression_factor)

    out = resample(data, downsample(data.Time, compression_factor));
end

