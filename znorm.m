function [ts] = znorm(ts)
    if (std(ts) ~= 0)
        ts = (ts - mean(ts)) ./ std(ts);
    else
        ts = (ts - ts(1));
    end
end

