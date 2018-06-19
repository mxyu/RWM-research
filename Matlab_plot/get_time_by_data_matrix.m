function [time, mat] = get_time_by_data_matrix( expt_time_data )
    %GET_TIME_BY_DATA_MATRIX argument
    time = unique(expt_time_data(:,1));
    mat = [];
    for i=1:length(time)
        data = expt_time_data(expt_time_data(:,1)==time(i),2); % data for a given time point
        mat(i,:) = data;
    end

end

