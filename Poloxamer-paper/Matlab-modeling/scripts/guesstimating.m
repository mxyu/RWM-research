filename = '/Users/michelle/Documents/2018_CUMC/Research/data_PBSandPOLOXAMER.xlsx';
[data_num,txt,raw] = xlsread(filename,2);
time = unique(data_num(:,2));

study = combine_expt_data(['H'], txt, data_num);
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

avg = [];
for i=1:size(study.data{1},1)
    avg(i) = mean(study.data{1}(i,:));
end

% plot(study.time, avg, 'b*')
x = study.time;
y = avg; % study.data{1}(:,6);
logx = log10(x(2:end));
logy = log10(y(2:end));

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
a = 0.0002; % in meters (200 microns)
k = 0; %thermal conductivity http://www.mhtl.uwaterloo.ca/pdf_papers/mhtl77-2.pdf
R_ss = 1/(4*k*a);
%D_arr = [10^-6];
D = 10^-6
factor = 10^3;
   
% for i=1:length(D_arr)
%     D = D_arr(i);
%     k = factor*D;
%     f = piecewise(D*t/a^2 < 0.6,... 
%         R_ss*(8/pi * (sqrt(D*t/a^2) - (D*t/a^2)/pi + (D*t/a^2)^2/(8*pi) + (D*t/a^2)^3/(32*pi) + 15*(D*t/a^2)^4/(512*pi))),...
%         D*t/a^2 >= 0.6,...
%         R_ss*(32/(3*pi^2) - (1/(pi^(3/2)*sqrt(D*t/a^2)))* (1 - 1/(3*(4*(D*t/a^2))) + 1/(6*(4*(D*t/a^2)^2)) - 1/(12*(4*(D*t/a^2)^3)))));
%     p(i) = fplot(f);
%     hold on;
% end

k = factor*D;
% f = piecewise(D*t/a^2 < 0.6,... 
%     R_ss*(8/pi * (sqrt(D*t/a^2) - (D*t/a^2)/pi + (D*t/a^2)^2/(8*pi) + (D*t/a^2)^3/(32*pi) + 15*(D*t/a^2)^4/(512*pi))),...
%     D*t/a^2 >= 0.6,...
%     R_ss*(32/(3*pi^2) - (1/(pi^(3/2)*sqrt(D*t/a^2)))* (1 - 1/(3*(4*(D*t/a^2))) + 1/(6*(4*(D*t/a^2)^2)) - 1/(12*(4*(D*t/a^2)^3)))));
R_ss = 1/(4*k*a);
f = symfun(R_ss*(8/pi * (sqrt(D*t/a^2) - (D*t/a^2)/pi + (D*t/a^2)^2/(8*pi) + (D*t/a^2)^3/(32*pi) + 15*(D*t/a^2)^4/(512*pi))), t);
f2 = symfun(R_ss*(32/(3*pi^2) - (1/(pi^(3/2)*sqrt(D*t/a^2)))* (1 - 1/(3*(4*(D*t/a^2))) + 1/(6*(4*(D*t/a^2)^2)) - 1/(12*(4*(D*t/a^2)^3)))), t);

p(1) = fplot(f);
hold on;
p(2) = fplot(f2);
hold on;

%legend(p, ["data polynomial" "10^-5" "10^-6" "10^-7"]);

p(3) = fplot(Rt_poly, 'Linewidth', 2);
%axis([0 170 0 3750])
ylabel("Resistance (min)", 'Fontsize', 20)
xlabel("Time (min)", 'Fontsize', 20)
set(gca,'fontsize',20)
hold on;

axis([0 165 0 2*10^6])

factor


