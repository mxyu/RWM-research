%% RUN ME FIRST: get data struct, polynomial model
filename = '/Users/michelle/Documents/2018_CUMC/Research/data_PBSandPOLOXAMER.xlsx';
[data_num,txt,raw] = xlsread(filename,2);
time = unique(data_num(:,2));

data_case = 'G';
% H = 1 large hole (200um diameter) in 18% Poloxamer
% G = 4 small holes (100um diameter) in 18% Poloxamer

study = combine_expt_data([data_case], txt, data_num);
% 
% p = [];
% for i = 1:size(study.data{1},2)
%     p(i) = plot(study.time, study.data{1}(:,i), ...
%         'LineWidth', 1);
%     hold on;
% end
% 
% expt_numbers = {};
% for i = 1:size(study.data{1},2)
%     expt_numbers{i} = int2str(i);
% end

% legend(p, expt_numbers)
% 
% avg = [];
% for i=1:size(study.data{1},1)
%     avg(i) = mean(study.data{1}(i,:));
% end
% 
% % plot(study.time, avg, 'b*')
% x = study.time;
% y = avg; % study.data{1}(:,6);
% logx = log10(x(2:end));
% logy = log10(y(2:end));

%% Polynomial model

b1 = 0.003808;
b2 = 0.58974;
b3 = 0.0037554;
syms x
% plot(study.time, avg, 'b*')
polyf = symfun(b1*x^b2 + b3, x);
% fplot(polyf, 'Linewidth', 2)
polydfdx = diff(polyf, x);
%% Plot current over time
% fplot(polydfdx, 'Linewidth', 2)
% axis([0 170 0 0.002])
% ylabel("Flow (concentration/mm^2/min)", 'Fontsize', 20)
% xlabel("Time (min)", 'Fontsize', 20)
% set(gca,'fontsize',20)
% Plot resistance over time
Rt_poly = 1/polydfdx;
fplot(Rt_poly, 'Linewidth', 2)
%axis([0 170 0 3750])
ylabel("Resistance (min)", 'Fontsize', 20)
xlabel("Time (min)", 'Fontsize', 20)
set(gca,'fontsize',20)
hold on;

%% Get diffusion constant by eye

% Model: Polynomial
p = [];


% Plot models with different magnitudes of diffusion constant
syms t
a = 0.0001; % in meters (diameter 200 microns)
D = 10^-9.5;
factor = 10^9;
k = factor*D;
R_ss = 1/(4*k*a);
t_bound = 0.6*a^2/D;
f = symfun(R_ss*(8/pi * (sqrt(D*t/(pi*a^2)) - (D*t/a^2)/pi + (D*t/a^2)^2/(8*pi) + (D*t/a^2)^3/(32*pi) + 15*(D*t/a^2)^4/(512*pi))), t);
f2 = symfun(R_ss*(32/(3*pi^2) - (2/(pi^(3/2)*sqrt(D*t/a^2)))* (1 - 1/(3*(4*(D*t/a^2))) + 1/(6*(4*(D*t/a^2))^2) - 1/(12*(4*(D*t/a^2))^3))), t);

% Create linspace for time to generate piecewise plot
time_min = 0;
time_max = 165;
time_points = 10^3;
time = linspace(time_min,time_max,time_points);
for i = 1:size(time,2)
    if time(i) < t_bound
        f_pw(i) = double(f(time(i))); % convert sym to double for pw array
    else
        f_pw(i) = double(f2(time(i)));
    end
end

p(1) = plot(time, f_pw, 'Linewidth', 2);
hold on;
p(2) = fplot(Rt_poly, 'Linewidth', 2);
hold on;

%axis([0 170 0 3750])
ylabel("Resistance (min)", 'Fontsize', 20)
xlabel("Time (min)", 'Fontsize', 20)
set(gca,'fontsize',20)

legend(p, ["f pw" "Rt poly"]);

%% Transform conductance model function from R(t) to f(t)
% Then we can fit it with the raw data and get a goodness of fit for the
% model
% g(t) = 1/R(t)
% int(g(t), t) = h(t)
% h(t) is in terms of concentration/mm^2


h_pw = [];
%time_test = linspace(0,165,10);

% Trapezoidal integration
f_pw_recip = 1./f_pw;
f_pw_recip_short = f_pw_recip(2:end);
time_short = time(2:end);
h_pw_test2 = cumtrapz(time_short,f_pw_recip_short);

time_interval = (time_max-time_min)/time_points;
for i = 1:size(time,2)
    h_pw_test(i) = (1/f_pw(i))*time_interval;
end

p(3) = plot(time_short, h_pw_test2, 'Linewidth', 2);
hold on;
%{
for i = 1:size(time,2)-1
    if time(i) < t_bound
        h_pw(i) = vpaintegral(1/f, t, [time(i) time(i+1)]); % convert sym to double for pw array
    else
        h_pw(i) = vpaintegral(1/f2, t, [time(i) time(i+1)]);
    end
end
%}

% p(1) = plot(time(1:end-1), h_pw, '--', 'Linewidth', 2);
% hold on;

% Plot raw study data
offset = 3;
for i = 1:size(study.data{1},2)
    p(i+offset) = plot(study.time, study.data{1}(:,i), ...
        'LineWidth', 1);
    hold on;
end

expt_numbers = {};
for i = 1:size(study.data{1},2)
    expt_numbers{i} = int2str(i);
end



%axis([0 170 0 3750])
ylabel("Concentration/mm^2", 'Fontsize', 20)
xlabel("Time (min)", 'Fontsize', 20)
set(gca,'fontsize',20)

legend(p(3:end), ["model" expt_numbers])


%% Resistance of 0 holes

[data_num,txt,raw] = xlsread(filename,2);
time = unique(data_num(:,2));

study = combine_expt_data(['F'], txt, data_num);
% Just used arbitrary area of 10


%% RUN ME SECOND: Total resistance accounting for rest of membrane (constant resistance)

area = 10;

% Polynomial fit
b1 = 0.003808;
b2 = 0.58974;
b3 = 0.0037554;
syms x
polyf = symfun(b1*x^b2 + b3, x);
%polydfdx = diff(polyf, x); 
% Convert to linspace
time_min = 0;
time_max = 165;
time_points = 10^3;
time = linspace(time_min,time_max,time_points);
time_interval = (time_max-time_min)/time_points;
ls_poly = [];
for i = 1:size(time,2)
    ls_poly(i) = double(polyf(time(i))); % convert sym to double for pw array
end

ls_poly2 = ls_poly(2:end);
ls_poly1 = ls_poly(1:end-1);
ls_polydiff = ls_poly2-ls_poly1;
ls_polyderiv = ls_polydiff/(time_interval)*10;

Rt_polyderiv = 1./ls_polyderiv;
% plot(time(1:end-1), Rt_polyderiv)


%Rt_poly = 1/polydfdx;
%% Play with coefficients of model terms

study = combine_expt_data(['H'], txt, data_num);

syms t
a = 0.0001; % in meters (diameter 200 microns)
D = 10^-10; %-15, -10
factor = 10^9.5;  % 13.1, 9.5?
k = factor*D;
R_ss = 1/(4*k*a);
t_bound = 0.6*a^2/D;
constant = 0.25;
first=0.25;
second=1.25;
third=1;
fourth=1;
fifth=1;
f = symfun(R_ss*constant*(first*8/pi * (sqrt(D*t/(pi*a^2)) - second*(D*t/a^2)/pi + third*(D*t/a^2)^2/(8*pi) + fourth*(D*t/a^2)^3/(32*pi) + fifth*15*(D*t/a^2)^4/(512*pi))), t);
f2 = symfun(R_ss*constant*(32/(3*pi^2) - (2/(pi^(3/2)*sqrt(D*t/a^2)))* (1 - 1/(3*(4*(D*t/a^2))) + 1/(6*(4*(D*t/a^2))^2) - 1/(12*(4*(D*t/a^2))^3))), t);
f_pw = [];
% Create linspace for time to generate piecewise plot
time_min = 0;
time_max = 165;
time_points = 10^3;
time = linspace(time_min,time_max,time_points);
for i = 1:size(time,2)
    if time(i) < t_bound
        f_pw(i) = double(f(time(i))); % convert sym to double for pw array
    else
        f_pw(i) = double(f2(time(i)));
    end
end


R_membrane = 335.0; % using arbitrary normalization area of 10mm^2 
R_total = 1./((1/R_membrane) + (1./f_pw));

% Plot resistance model with data fit
%
f1 = figure();
p=[];
p(1) = plot(time, R_total,'Linewidth', 2);
hold on;
p(2) = plot(time(1:end-1), Rt_polyderiv,'Linewidth', 2);
hold on;
legend(p, ["model" "data fit"]);
%

% Trapezoidal integration
R_total_recip = 1./R_total;
R_total_recip_short = R_total_recip(2:end);
time_short = time(2:end);
R_total_integral = cumtrapz(time_short,R_total_recip_short);
R_integral_normalized = R_total_integral/10;

% Test trapezoidal integration on Rtpolyderiv
Rt_polyderiv_recip = 1./Rt_polyderiv;
Rt_integral = cumtrapz(time_short,Rt_polyderiv_recip);
Rt_integral_normalized = Rt_integral/10;

% Plot final model with raw data
f2 = figure;
p=[];
p(1) = plot(time_short, R_integral_normalized, 'Linewidth', 2);
hold on;
p(2) = plot(time_short, Rt_integral_normalized, 'Linewidth', 2);
hold on;
legend(p, ["model" "data fit"]);

offset = 3;
for i = 1:size(study.data{1},2)
    p(i+offset) = plot(study.time, study.data{1}(:,i), ...
        'LineWidth', 0.25);
    hold on;
end

%store model in new variable name so we can compare with 4 hole model

model_1hole = R_integral_normalized;


%% Modeling 4 holes data 

study = combine_expt_data(['G'], txt, data_num);

syms t
a = 0.0001; % in meters (diameter 100 microns, a is radius)
D = 10^-10; %-15, -10
factor = 10^9.5;  % 13.1, 9.5?
k = factor*D;
R_ss = 1/(4*k*a);
t_bound = 0.6*a^2/D;
constant = 0.25;
first=0.25;
second=1.25;
third=1;
fourth=1;
fifth=1;
f = symfun(R_ss*constant*(first*8/pi * (sqrt(D*t/(pi*a^2)) - second*(D*t/a^2)/pi + third*(D*t/a^2)^2/(8*pi) + fourth*(D*t/a^2)^3/(32*pi) + fifth*15*(D*t/a^2)^4/(512*pi))), t);
f2 = symfun(R_ss*constant*(32/(3*pi^2) - (2/(pi^(3/2)*sqrt(D*t/a^2)))* (1 - 1/(3*(4*(D*t/a^2))) + 1/(6*(4*(D*t/a^2))^2) - 1/(12*(4*(D*t/a^2))^3))), t);
f_pw = [];
% Create linspace for time to generate piecewise plot
time_min = 0;
time_max = 165;
time_points = 10^3;
time = linspace(time_min,time_max,time_points);
for i = 1:size(time,2)
    if time(i) < t_bound
        f_pw(i) = double(f(time(i))); % convert sym to double for pw array
    else
        f_pw(i) = double(f2(time(i)));
    end
end


R_membrane = 335.0; % using arbitrary normalization area of 10mm^2 
R_total = 1./((1/R_membrane) + 4*(1./f_pw)); % Treating holes as resistors in parallel

% Plot resistance model with data fit
%
f1 = figure();
p=[];
p(1) = plot(time, R_total,'Linewidth', 2);
hold on;
p(2) = plot(time(1:end-1), Rt_polyderiv,'Linewidth', 2);
hold on;
legend(p, ["model" "data fit"]);
%

% Trapezoidal integration
R_total_recip = 1./R_total;
R_total_recip_short = R_total_recip(2:end);
time_short = time(2:end);
R_total_integral = cumtrapz(time_short,R_total_recip_short);
R_integral_normalized = R_total_integral/10;

% Test trapezoidal integration on Rtpolyderiv
Rt_polyderiv_recip = 1./Rt_polyderiv;
Rt_integral = cumtrapz(time_short,Rt_polyderiv_recip);
Rt_integral_normalized = Rt_integral/10;

% Plot final model with raw data
f2 = figure;
p=[];
p(1) = plot(time_short, R_integral_normalized, 'Linewidth', 2);
hold on;
p(2) = plot(time_short, Rt_integral_normalized, 'Linewidth', 2);
hold on;
legend(p, ["model" "data fit"]);

offset = 3;
for i = 1:size(study.data{1},2)
    p(i+offset) = plot(study.time, study.data{1}(:,i), ...
        'LineWidth', 0.25);
    hold on;
end

model_4holes = R_integral_normalized;

%% Plot 1 hole vs 4 hole models

p = [];
p(1) = plot(time_short, model_1hole);
hold on;
p(2) = plot(time_short, model_4holes);
hold on;
legend(p,["1 hole model" "4 hole model"]);

%% Plot models with different radius and different hole number

% Different hole radius
r = [50 50 50];

% Hole numbers
h = [4 8 16];

% Membrane resistance
m = [335 335 335];

mat = [r' h' m'];


models = {};
for i = 1:size(mat,1)
    models{i} = generate_model(mat(i,:));
end

p = [];
for j = 1:size(mat,1)
    p(j) = plot(time_short, models{j});
    hold on;
end

labels = {};
for k = 1:size(mat,1)
    labels{k} = "diameter: " + mat(k,1)*2 + ", # holes: " + mat(k,2) + ", R_mem: " + mat(k,3);
end

%legend(p, labels);

