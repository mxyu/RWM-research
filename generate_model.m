function model = generate_model( vars )
%GENERATE_MODEL Summary of this function goes here
%   radius = radius of hole, given in args in microns. Convert by *10^-6
%   R_membrane = resistance of membrane (from experimental measures)

radius = vars(1);
num_holes = vars(2);
R_membrane = vars(3);


syms t
a = radius*10^-6; % in meters (diameter 100 microns, a is radius)
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


%R_membrane = 335.0; % using arbitrary normalization area of 10mm^2 
R_total = 1./((1/R_membrane) + num_holes*(1./f_pw)); % Treating holes as resistors in parallel

% Trapezoidal integration
R_total_recip = 1./R_total;
R_total_recip_short = R_total_recip(2:end);
time_short = time(2:end);
R_total_integral = cumtrapz(time_short,R_total_recip_short);
R_integral_normalized = R_total_integral/10;


model = R_integral_normalized;

end

