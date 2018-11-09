% Modelling Jaundice
clc; 
clear all;
C1 = 0.25 ; % [mg/dl/hr] - hourly creation of bilirubin by the breakdown of red blood cells
a = [0 0 1; 0 1 1; 1 1 0; 1 0 1]; % Constant of proportionality
syms x1(t) x2(t) x3(t)

eqn1 = diff(x1) == -(a(2,1) + a(3,1) + a(4,1))*x1 + a(1,2)*x2 + a(1,3)*x3 + C1
eqn2 = diff(x2) == a(2,1)*x1 - (a(1,2) + a(4,2))*x2
eqn3 = diff(x3) == a(3,1)*x1 - a(1,3) + a(4,3)*x3

c1 = x1(0) == 1; % initial bilirubin concentrations in the body surface
c2 = x2(0) == 1; % initial bilirubin concentrations in the blood 
c3 = x3(0) == 1; % initial bilirubin concentrations in the liver

[x1Sol(t) x2Sol(t) x3Sol(t)] = dsolve(eqn1, eqn2, eqn3, c1, c2, c3)

fplot(x1Sol)
hold on
fplot(x2Sol)
hold on
fplot(x3Sol)
grid on
legend('x1Sol','x2Sol','x3Sol','Location','best')