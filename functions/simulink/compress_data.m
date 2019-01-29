function [new_data] = compress_data(data, compression_rate)
%COMPRESS_DATA Summary of this function goes here
%   Detailed explanation goes here
    
    if nargin < 2
        compression_rate = 1e3;
    end
    
    new_data = data(1:compression_rate:end);
end

