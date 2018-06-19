function [time, xmean, yerr, ystd] = get_mean_err_by_time(expt_time_data)
    time = unique(expt_time_data(:,1));
    xmean = [];
    yerr = [];
    ystd = [];
    for i=1:length(time)
        data = expt_time_data(expt_time_data(:,1)==time(i),2);
        xmean(end+1) = mean(data);
        ystd(end+1) = std(data);
        se = std(data) / sqrt(length(data));
        yerr(end+1) = 2*se;
    end
end