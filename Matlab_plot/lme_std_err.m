function CI = lme_std_err( lme )
%LME_STD_ERR Summary of this function goes here
%   Detailed explanation goes here
    NUM_OBSERVATIONS = 16;
    CI.x = [0:1:165];
    CI.std = [];
    for i = 1:length(CI.x)
        CI.std(end+1) = sqrt(CI.x(i)^2 * lme.s_err^2 + lme.i_err^2) * 2 / sqrt(NUM_OBSERVATIONS);
    end
    CI.upper = lme.slope * CI.x + lme.intercept + CI.std;
    CI.lower = lme.slope * CI.x + lme.intercept - CI.std;

end

