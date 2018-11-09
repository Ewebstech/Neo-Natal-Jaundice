% Calculating Liver Performance of a Human Patient
clc;
clear all;
close all;

hrs = 120;
P0 = 0.05;              % Ability to excrete bilirubin
%matrate = 0.;         % Maturity Rate of Liver
P0 = 0.05;              % Ability to excrete bilirubin
r_m = 0.7;    % Maturity Rate Vector
time = 1:1:hrs;

for t = 1:length(time)
    P(t) = P0/(P0 + (1-P0)*exp(-r_m*t))
 
end
Percentage = P .* 100; %
p1 = plot(time,Percentage)
p1.LineWidth = 2;

hold on

hrs = 120;
matrate = 0.01; % Maturity Rate of Liver
P0 = 0.05; % Ability to excrete bilirubin
r_m = 0.09; % Maturity Rate Vector
time = 1:1:hrs;

for t = 1:length(time)
    P(t) = P0/(P0 + (1-P0)*exp(-r_m*t))
 
end
Percentage = P .* 100; %
p2 = plot(time,Percentage,':')
p2.LineWidth = 1;
p2.Marker = '*';

xlabel('Time (Hrs)'), ylabel('Percentage of Liver Performance')
title('Liver Performance Graph'), grid on,
legend('4% Mat. Rate r_m = 0.7', '%0.1 Mat. Rate r_m = 0.09')

