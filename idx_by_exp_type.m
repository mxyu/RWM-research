function idx_arr = idx_by_exp_type(data_txt,expType)
    tmp = strfind(data_txt(:,1),expType);
    idx_arr_tmp = zeros([length(tmp)-1,1]);
    for i = 1:length(idx_arr_tmp)
        if int8(tmp{i+1}) == 1
            idx_arr_tmp(i) = 1;
        end
    end
    idx_arr = find(idx_arr_tmp == 1);
end