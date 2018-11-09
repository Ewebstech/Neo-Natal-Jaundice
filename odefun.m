function dxdt = odefun(t,x,hrs)
global LiverP n_time C1 CT_1 SnMP Blood_Transfusion Phototherapy
dxdt = zeros(3,1);

% Computational Experiment has shown that vector 'a' should not be > 1
 xmin = 1;
 xmax = 1;
 n = 3; m = 4;

% Unitary Value Converges

a = xmin+rand(4,3)*(xmax-xmin); % Constant of Proportionality
%a = [0, 45.123, 40.456; 90.789, 0, 0; 20.007, 0, 0; 28.099, 31.055, 29.567]
%disp('Computing and Calculating Differential Systems...');
% Calculate a13 and a43 Liver Performace Data
if (n_time < hrs)
    n_time = n_time + 1;
end

a13 = 1/sqrt(LiverP(n_time))*a(1,3);
a43 = a(4,3)*LiverP(n_time);


% Treatment by Blood Transfusion (15-20 mg/dl)
C0 = 1.0; % mg/dl Bilirubin Conc of New Blood
r_blood = 0.002; % l/kg/hr Transfusion Rate
V = 0.000109; %micro-mL/kg Volume of Blood in the human body

% Treatment by Phototherapy (10-15md/dl)
k12 = 6; % micro-watts/cm^2/nm Intensity of Bililight for compartment 1 and 2
k42 = 6; % micro-watts/cm^2/nm; % Intensity of Bililight for compartment 4 and 2

if(SnMP == 1)
    C = CT_1;
else
    C = C1;
end

% System of Differential Equations
dxdt(1) = -(a(2,1) + a(3,1) + a(4,1))*x(1) + a(1,2)*x(2) + x(3)*a13 + C(n_time) + (r_blood/V)*(C0 - x(1))*Blood_Transfusion;
dxdt(2) = a(2,1)*x(1) - (a(1,2) + a(4,2))*x(2) - (a(1,2)*k12 + a(4,2)*k42)*x(2)*Phototherapy;
dxdt(3) = a(3,1)*x(1) - (a13 + a43)*x(3);


