function tmp = generate_lme_plots(study, colors)
%GENERATE_LME_PLOTS Summary of this function goes here
%   Detailed explanation goes here

    % Plot LME in panels
    f = figure('Name',study.title);
    for i = 1:length(study.data)
        subplot(1,length(study.data), i);
        plot(study.time, study.data{i}, ...
            'LineWidth', 1, ...
            'color', colors.GRAY);
        title(study.expt_names{i})
        axis(study.axis)
        hold on
        plot(study.time, study.lmes{i}, ...
            'LineWidth', 2, ...
            'color', colors.BLUE);
        xlabel(study.x_axis_label,'FontSize',12)
        if i == 1
            ylabel(study.y_axis_label,'FontSize',12)
        end
        set(gca,'box','off')
        set(gca, 'TickLength',[0 0])
    end
    
    
    f2 = figure('Name',study.title);
    p = [];
    for i = 1:length(study.data)
        plot(study.time, study.data{i}, ...
            'LineWidth', 1, ...
            'color', colors.LIGHT{i});
        hold on
    end
    for i = 1:length(study.data)
        p(i) = plot(study.time, study.lmes{i}, ...
            'LineWidth', 2, ...
            'color', colors.PRIMARY{i});
        hold on
    end
    legend(p, study.expt_names,'Location','northwest')
    xlabel(study.x_axis_label,'FontSize',12)
    ylabel(study.y_axis_label,'FontSize',12)
    axis(study.axis)
    title(study.title)
    set(gca,'box','off')
    set(gca, 'TickLength',[0 0])
    hold on
end

