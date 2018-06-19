function study = combine_expt_data(expt_arr, txt, data_num)
    OUTCOME_COL = 7;    
    study = {};
    study.data = {};
    study.time = [];    
    for i=1:length(expt_arr)
        [t, d] = get_time_by_data_matrix(data_num(idx_by_exp_type(txt,expt_arr(i)),[2 OUTCOME_COL]));
        study.time = t;
        study.data{i} = d;
    end
    
end

