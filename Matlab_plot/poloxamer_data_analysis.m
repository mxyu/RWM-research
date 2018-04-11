filename = '/Users/michelle/Documents/2018_CUMC/Research/data_PBSandPOLOXAMER.xlsx';
[data_num,txt,raw] = xlsread(filename,2);

OUTCOME_COL = 7;
GRAY = [230, 230, 230]/255;
RED = [96.1, 0, 16.9]/100;
GREEN = [24.7, 83.1, 0]/100;
BLUE = [4.3, 47.5, 78]/100;
ORANGE = [100, 56.5, 0]/100;
PURPLE = [26.7, 8.6, 80.8]/100;
YELLOW = [100, 83.5, 0]/100;

LRED = [100, 83.9, 86.7]/100;
LGREEN = [88.6, 100, 83.9]/100;
LBLUE = [85.1, 93.7, 100]/100;
LORANGE = [100, 92.9, 83.9]/100;
LPURPLE = [89.4, 85.9, 100]/100;
LYELLOW = [100, 83.5, 0]/100;

COLORS = {};
COLORS.GRAY = GRAY;
COLORS.BLUE = BLUE;
COLORS.PRIMARY = {RED; BLUE; GREEN; ORANGE; YELLOW; PURPLE};
COLORS.LIGHT = {LRED; LBLUE; LGREEN; LORANGE; LYELLOW; LPURPLE};

time = unique(data_num(:,2));


%% Plot Study 1: A-D, Comparing different # holes, same hole diam, PBS
% LME
lme_ABCD_intercepts = [ -0.005721924, 0.020320917+ -0.005721924, 0.015407894+ -0.005721924, 0.017743872+ -0.005721924];
lme_ABCD_slopes = [0.003176679, 0.000362071+0.003176679, 0.000865643+0.003176679, 0.001505025+0.003176679];

A_lme = lme_ABCD_intercepts(1) + time * lme_ABCD_slopes(1);
B_lme = lme_ABCD_intercepts(2) + time * lme_ABCD_slopes(2);
C_lme = lme_ABCD_intercepts(3) + time * lme_ABCD_slopes(3);
D_lme = lme_ABCD_intercepts(4) + time * lme_ABCD_slopes(4);

study1 = combine_expt_data(['A','B','C','D'], txt, data_num);
study1.title = 'Concentration/mm^2 over Time for 0, 1, 2, and 4 Holes in PBS';
study1.lmes = {A_lme ; B_lme ; C_lme; D_lme};
study1.expt_names = {'0 holes'; '1 hole'; '2 holes'; '4 holes'};
study1.axis = [0 170 0 1.15];
study1.x_axis_label = {'Time (min)'};
study1.y_axis_label = {'Concentration/mm^2'};

generate_lme_plots(study1, COLORS);

%% Plot Study 2: A,D,E, Comparing different # holes, same hole AREA, PBS
A_lme = -0.005721924 + time * 0.003176679;
D_lme = -0.005721924 + 0.017743872 + time * ( 0.003176679 + 0.001505025);  
E_lme = -0.005721924 - 0.011300950 + time * ( 0.003176679 + 0.001567684);

study2 = combine_expt_data(['A','E','D'], txt, data_num);
study2.title = 'Concentration/mm^2 over Time for 0, 1, and 4 Holes in PBS';
study2.lmes = {A_lme ; D_lme; E_lme};
study2.expt_names = {'0 holes'; '1 hole, 200 um diameter'; '4 holes, 100 um diameter'};
study2.axis = [0 170 0 1];
study2.x_axis_label = {'Time (min)'};
study2.y_axis_label = {'Concentration/mm^2'};

generate_lme_plots(study2, COLORS);


%% Plot study 3 - Compare all poloxamer data points (F, G, H)
F_lme = 0.008449670 + time * 0.000298481;
H_lme = 0.008449670 + 0.007890833 + time * (0.000298481 + 0.000125602); 
G_lme = 0.008449670 + 0.003805651 + time * (0.000298481 + 0.000023919); 

study3 = combine_expt_data(['F','G','H'], txt, data_num);
study3.title = 'Concentration/mm^2 over Time for 0, 1, and 4 Holes in 18% Poloxamer';
study3.lmes = {F_lme ; G_lme; H_lme};
study3.expt_names = {'0 holes'; '1 hole, 200 um diameter'; '4 holes, 100 um diameter'};
study3.axis = [0 170 0 0.12];
study3.x_axis_label = {'Time (min)'};
study3.y_axis_label = {'Concentration/mm^2'};

generate_lme_plots(study3, COLORS);


%% Plot study 4 - Compare PBS same area expt to Poloxamer 18% same area expt (ADE vs FGH)
lme_A = -0.005721924 + time * 0.003176679;
lme_D = (-0.005721924 + 0.017743872) + time * (0.003176679 + 0.001505025);
lme_E = (-0.005721924 + -0.011300950) + time * (0.003176679 + 0.001567684);
lme_F = (-0.005721924 + 0.014171594) + time * (0.003176679 + -0.002878198);
lme_G = (-0.005721924 + 0.017977245) + time * (0.003176679 + -0.002854279);
lme_H = (-0.005721924 + 0.022062427) + time * (0.003176679 + -0.002752596);

study = combine_expt_data(['A','E','D','F','H','G'], txt, data_num);
study.title = 'Concentration/mm^2 over Time for 0, 1, and 4 Holes in PBS and 18% Poloxamer';
study.lmes = {A_lme; E_lme; D_lme; F_lme ; H_lme; G_lme};
study.expt_names = {'PBS: 0 holes'; 'PBS: 1 hole, 200 um diameter'; 'PBS: 4 holes, 100 um diameter'; '18% Poloxamer: 0 holes'; '18% Poloxamer: 1 hole, 200 um diameter'; '18% Poloxamer: 4 holes, 100 um diameter'};
study.axis = [0 170 0 0.8];
study.x_axis_label = {'Time (min)'};
study.y_axis_label = {'Concentration/mm^2'};

f2 = figure('Name',study.title);
p = [];
for i = 1:3
    p(i) = plot(study.time, study.lmes{i}, ...
        'LineWidth', 1, ...
        'color', COLORS.PRIMARY{i});
    hold on
end
for i = 4:6
    p(i) = plot(study.time, study.lmes{i}, '--', ...
        'LineWidth', 1, ...
        'color', COLORS.PRIMARY{i-3});
    hold on
end
legend(p, study.expt_names,'Location','northwest')
xlabel(study.x_axis_label,'FontSize',12)
ylabel(study.y_axis_label,'FontSize',12)
axis(study.axis)
title(study.title)
set(gca,'box','off')
set(gca, 'TickLength',[0 0])

%% Plot 0 holes, PBS vs 18% Poloxamer

lme_A = -0.005721924 + time * 0.003176679;
lme_F = (-0.005721924 + 0.014171594) + time * (0.003176679 + -0.002878198);

lme2 = {};
lme2.intercept = -0.005721924;
lme2.slope = 0.003176679;
lme_ref = {};
lme_ref.intercept = (-0.005721924 + 0.014171594);
lme_ref.slope = (0.003176679 + -0.002878198);

b1 = plot_comparison_line_mat(lme_ref, lme2, 60);
b2 = plot_comparison_line_mat(lme_ref, lme2, 120);

study = combine_expt_data(['A','F'], txt, data_num);
study.title = 'Concentration/mm^2 over Time for 0 Holes in PBS vs 18% Poloxamer';
study.lmes = {lme_A; lme_F };
study.expt_names = {'PBS: 0 holes'; '18% Poloxamer: 0 holes'};
study.axis = [0 170 0 0.15];
study.x_axis_label = {'Time (min)'};
study.y_axis_label = {'Concentration/mm^2'};

f2 = figure('Name',study.title);
p = [];
for i = 1:length(study.data)
    plot(study.time, study.data{i}, ...
        'LineWidth', 1, ...
        'color', COLORS.LIGHT{i});
    hold on
end
for i = 1:length(study.data)
    p(i) = plot(study.time, study.lmes{i}, ...
        'LineWidth', 2, ...
        'color', COLORS.PRIMARY{i});
    hold on
end

for i = 1:size(b1,3)
    plot(b1(1,:,i), b1(2,:,i), ':', 'color', 'black', 'LineWidth', 1.5)
    hold on
end
for i = 1:size(b2,3)
    plot(b2(1,:,i), b2(2,:,i), ':', 'color', 'black', 'LineWidth', 1.5)
    hold on
end

legend(p, study.expt_names,'Location','northwest')
xlabel(study.x_axis_label,'FontSize',12)
ylabel(study.y_axis_label,'FontSize',12)
axis(study.axis)
title(study.title)
set(gca,'box','off')
set(gca, 'TickLength',[0 0])


%% Plot 1 hole, PBS vs 18% Poloxamer

lme_E = -0.01702287 + time * 0.00474436;
lme_H = (-0.01702287 + 0.03336338) + time * ( 0.00474436 + -0.00432028);

lme2 = {};
lme2.intercept = -0.01702287;
lme2.slope = 0.00474436;
lme_ref = {};
lme_ref.intercept = -0.01702287 + 0.03336338;
lme_ref.slope = 0.00474436 + -0.00432028;

b1 = plot_comparison_line_mat(lme_ref, lme2, 60);
b2 = plot_comparison_line_mat(lme_ref, lme2, 120);

study = combine_expt_data(['E','H'], txt, data_num);
study.title = 'Concentration/mm^2 over Time for 1 Hole in PBS vs 18% Poloxamer';
study.lmes = {lme_E; lme_H };
study.expt_names = {'PBS: 1 hole'; '18% Poloxamer: 1 hole'};
study.axis = [0 170 0 0.2];
study.x_axis_label = {'Time (min)'};
study.y_axis_label = {'Concentration/mm^2'};

f2 = figure('Name',study.title);
p = [];
for i = 1:length(study.data)
    plot(study.time, study.data{i}, ...
        'LineWidth', 1, ...
        'color', COLORS.LIGHT{i});
    hold on
end
for i = 1:length(study.data)
    p(i) = plot(study.time, study.lmes{i}, ...
        'LineWidth', 2, ...
        'color', COLORS.PRIMARY{i});
    hold on
end

for i = 1:size(b1,3)
    plot(b1(1,:,i), b1(2,:,i), ':', 'color', 'black', 'LineWidth', 1.5)
    hold on
end
for i = 1:size(b2,3)
    plot(b2(1,:,i), b2(2,:,i), ':', 'color', 'black', 'LineWidth', 1.5)
    hold on
end

legend(p, study.expt_names,'Location','northwest')
xlabel(study.x_axis_label,'FontSize',12)
ylabel(study.y_axis_label,'FontSize',12)
axis(study.axis)
title(study.title)
set(gca,'box','off')
set(gca, 'TickLength',[0 0])


%% Scratch work

% % Plot linear regression
% %Plot A
% pA = plot(time, A_lme, ...
%     'LineWidth', 2, ...
%     'color', RED);
% hold on
% % Plot B
% pB = plot(time, B_lme, ...
%     'LineWidth', 2, ...
%     'color', ORANGE);
% hold on
% % Plot C
% pC = plot(time, C_lme, ...
%     'LineWidth', 2, ...
%     'color', GREEN);
% hold on
% % Plot D
% pD = plot(time, D_lme, ...
%     'LineWidth', 2, ...
%     'color', BLUE);
% legend([pA pB pC pD], '0 holes', '1 hole', '2 holes', '4 holes','Location','northwest')
% xlabel('Time (minutes)','FontSize',12)
% ylabel('Concentration/mm^2','FontSize',12)
% title('Concentration/mm^2 over Time for 0, 1, 2, and 4 Holes in PBS')
% %plot(T.Concentration, T.numHoles, 'ro')
% %xlabel('Concentration')
% %ylabel('Num Holes')



% Compare 0 holes, 1 hole, 2 hole, 4 holes (all 100um) vs Concentration
% RhoB
% Compare PBS, Poloxamer vs Concentration RhoB 
% Compare 2 100um holes, vs 1 200um hole vs Concentration RhoB

% % Plot A-D
% % mean, standard error * 2
% 
% A_time_by_conc = data_num(idx_by_exp_type(txt,'A'),[2 OUTCOME_COL]);
% %  mean and std err at each time point
% [time, A_means, A_err, A_std] = get_mean_err_by_time(A_time_by_conc);
% 
% % Get mean, standard error, and standard deviation of experiments A-D
% [time, B_means, B_err, B_std] = get_mean_err_by_time(data_num(idx_by_exp_type(txt,'B'),[2 OUTCOME_COL]));
% [time, C_means, C_err, C_std] = get_mean_err_by_time(data_num(idx_by_exp_type(txt,'C'),[2 OUTCOME_COL]));
% [time, D_means, D_err, D_std] = get_mean_err_by_time(data_num(idx_by_exp_type(txt,'D'),[2 OUTCOME_COL]));
% [time, E_means, E_err, E_std] = get_mean_err_by_time(data_num(idx_by_exp_type(txt,'E'),[2 OUTCOME_COL]));

% f = figure('Name','Experiments in PBS');
% % plot(time, A_means,...
% %     'LineWidth',2,...
% %     'color', [179, 224, 255]/255);
% hold on
% errorbar(time, A_means, A_err,...
%       'LineWidth',1,...
%       'color', [179, 224, 255]/255);
% hold on;
% errorbar(time, B_means, B_err,...
%       'LineWidth',1,...
%       'color', [26, 163, 255]/255);
% hold on;
% errorbar(time, C_means, C_err,...
%       'LineWidth',1,...
%       'color', [0, 122, 204]/255);
% hold on;
% errorbar(time, D_means, D_err,...
%       'LineWidth',1,...
%       'color', [0, 77, 128]/255);
% %legend('Exp A', 'Exp B', 'Exp C', 'Exp D','Location','northwest')

