function box = plot_comparison_line_mat(lme_ref, lme2, time)
    %UNTITLED Summary of this function goes here
    %   Detailed explanation goes here
    lme_ref_at_time = lme_ref.intercept + (time) * lme_ref.slope;
    %Find time at which PBS concentration is equal to Polo concentration at 1hr
    %and 2hr
    lme2_at_time = (lme_ref_at_time - lme2.intercept) / lme2.slope;

    % Line coordinates
    box(:,:,1) = [time time; -1 lme_ref_at_time];
    box(:,:,2) = [lme2_at_time time; lme_ref_at_time lme_ref_at_time];
    box(:,:,3) = [lme2_at_time lme2_at_time; -1 lme_ref_at_time];

end

