% Modelling Jaundice
clc; 
clear all;
global SnMP C1 CT_1 C2 Phototherapy Blood_Transfusion; % Declaring Value of C as a Global Variable
hrs = 120; % Number of experimental hours
%% %Treatment Options
Blood_Transfusion = 0;
Phototherapy = 0;
SnMP = 0;

P0 = 0.05; % Ability to excrete bilirubin
r_m = 0.09; % mat rate of liver

bil_rate_no_treatment = 0.9; %mg/dl [mg/dl/hr] - hourly creation of bilirubin by the breakdown of red blood cells
bil_rate_treatment = 0.05; %mg/dl [mg/dl/hr] - New hourly creation of bilirubin by the breakdown of red blood cells

C1_rate = 0;
for j = 1:hrs
    C1_rate = C1_rate + bil_rate_no_treatment;
    C1(j) = C1_rate;
end

C2_rate = 0;
for r = 1:hrs
    C2_rate = C2_rate + bil_rate_treatment;
    C2(r) = C2_rate;
end

T_3 = 5; % mol/kg (6 micromol/kg) Drug Dosage
D = 6; % Hours of Drug Administration
for n = 1:hrs
    val = ((n - T_3)*pi)/(D);
    CT_1(n) = C2(n)  + (C1(n) - C2(n))*cos(val)*SnMP;
end
tspan = [0 8];          % Time Interval

% Find the percentage of normal Liver function performance
time = 1:1:hrs;

for i = 1:length(time)
    Perc(i) = P0/(P0 + (1-P0)*exp(-r_m*i));
end
P = Perc;

global LiverP n_time;
LiverP = P;
n_time = 0;
x0 = [1.85 0.28 5.40]; % initial bilirubin concentrations in the body surface
% Using Ode113 Matlab Function ODE Solve
opts = odeset('Reltol',1e-6,'AbsTol',1e-6,'Stats','on');
[t,x] = ode113(@(t,x) odefun(t,x,hrs), tspan, x0, opts);

x1_t = x(40:end,1); 
x2_t = x(40:end,2);
x3_t = x(40:end,3);
timeT = 40:1:length(x(:,2));
time_T = timeT(:);
%clc;
table(time_T,x1_t,x2_t,x3_t)
% Perform Plots
%p = plot(time_T,x(:,1),'k--',time_T,x(:,2),'r',time_T,x(:,3));
figure 
subplot(2,1,1);
p = plot(time_T,x1_t,'k:');
p.LineWidth = 2;
hold on
q = plot(time_T,x2_t,'r');
q.LineWidth = 2;
hold on
r = plot(time_T,x3_t,'g-.');
r.LineWidth = 2;

legend('Body Surface Conc.', 'Blood Conc.', 'Liver Conc.');
xlabel('Time(hrs)'), ylabel('Concentration mg/dl')
title('Bilirubin Blood Concentration Level'), grid on
hold on


%% Treatments
Blood_Transfusion = input('Switch on Treatment by Blood Transfusion: Enter Value >> ');
Phototherapy = input('Switch on Treatment by Phototherapy: Enter Value >> ');
SnMP = input('Switch on Treatment by Sn-MP Medication: Enter Value >> ');

%PT = x2_t(find(x2_t > 10));
XT1_value = x1_t(end);
XT2_value = x2_t(end);
XT3_value = x3_t(end);

Time_end = time_T(end);

if(Blood_Transfusion == 1)
    Treatment(2) = 1;
    x0 = [XT1_value XT2_value XT3_value]; % initial bilirubin concentrations
    % Using Ode113 Matlab Function ODE Solve
    opts = odeset('Reltol',1e-6,'AbsTol',1e-6,'Stats','on');
    [t,x] = ode113(@(t,x) odefun(t,x,hrs), tspan, x0, opts);

    x1_t = x(:,1);
    x2_t = x(:,2);
    x3_t = x(:,3);
    timeT = Time_end:1:Time_end+(length(x2_t)-1);
    time_T = timeT(:);
    
    table(time_T,x1_t,x2_t,x3_t)
    % Perform Plots
    p = plot(time_T,x1_t,'k:');
    p.LineWidth = 2;
    hold on
    q = plot(time_T,x2_t,'r');
    q.LineWidth = 2;
    hold on
    r = plot(time_T,x3_t,'g-.');
    r.LineWidth = 2;
    
end

if(Phototherapy == 1) 
    Treatment(1) = 1;
    x0 = [XT1_value XT2_value XT3_value]; % initial bilirubin concentrations
    % Using Ode113 Matlab Function ODE Solve
    opts = odeset('Reltol',1e-6,'AbsTol',1e-6,'Stats','on');
    [t,x] = ode113(@(t,x) odefun(t,x,hrs), tspan, x0, opts);

    x1_t = x(:,1);
    x2_t = x(:,2);
    x3_t = x(:,3);
    timeT = Time_end:1:Time_end+(length(x2_t)-1);
    time_T = timeT(:);
    
    table(time_T,x1_t,x2_t,x3_t)
    % Perform Plots
    p = plot(time_T,x1_t,'k:');
    p.LineWidth = 2;
    hold on
    q = plot(time_T,x2_t,'r');
    q.LineWidth = 2;
    hold on
    r = plot(time_T,x3_t,'g-.');
    r.LineWidth = 2;
   
end

if(SnMP == 1) 
    Treatment(3) = 1;
    x0 = [XT1_value XT2_value XT3_value]; % initial bilirubin concentrations
    % Using Ode113 Matlab Function ODE Solve
    opts = odeset('Reltol',1e-6,'AbsTol',1e-6,'Stats','on');
    [t,x] = ode113(@(t,x) odefun(t,x,hrs), tspan, x0, opts);

    x1_t = x(:,1);
    x2_t = x(:,2);
    x3_t = x(:,3);
    timeT = Time_end:1:Time_end+(length(x2_t)-1);
    time_T = timeT(:);
    clc;
    table(time_T,x1_t,x2_t,x3_t)
    % Perform Plots
    p = plot(time_T,x1_t,'k:');
    p.LineWidth = 2;
    hold on
    q = plot(time_T,x2_t,'r');
    q.LineWidth = 2;
    hold on
    r = plot(time_T,x3_t,'g-.');
    r.LineWidth = 1;
   
end
if(SnMP == 1) 
    subplot(2,1,2);
    time = 1:1:length(x2_t);
    time_t = time(:);
    rt = (x1_t + x2_t + x3_t)/(2^1);
    tr_plot = plot(time_T,rt,'k-.'); 
    xlabel('Time(hrs)'), ylabel('Treatment')
    hold on
end
if(Phototherapy == 1) 
    subplot(2,1,2);
    time = 1:1:length(x2_t);
    time_t = time(:);
    rt = (x1_t + x2_t + x3_t)/(2^2);
    tr_plot = plot(time_T,rt,'g--'); 
    xlabel('Time(hrs)'), ylabel('Treatment')
    hold on
end
if(Blood_Transfusion == 1) 
    subplot(2,1,2);
    time = 1:1:length(x2_t);
    time_t = time(:);
    rt = (x1_t + x2_t + x3_t)/(2^3);
    tr_plot = plot(time_T,rt,'r-*'); 
    xlabel('Time(hrs)'), ylabel('Treatment')
end

if(Phototherapy == 1 && SnMP == 0 && Blood_Transfusion == 0)
    title('Bilirubin Conc. Treatment with Phototherapy Alone');
    legend('Phototherapy');
elseif(SnMP == 1 && Phototherapy == 0 && Blood_Transfusion == 0)
    title('Bilirubin Conc. Treatment with Sn-MP Medication Alone');
    legend('Sn-MP');
elseif(Blood_Transfusion == 1 && SnMP == 0 && Phototherapy == 0)
    title('Bilirubin Conc. Treatment with Blood Transfusion Alone');
    legend('Blood Transfusion');
elseif(Blood_Transfusion == 1 && SnMP == 1 && Phototherapy == 0)
    title('Blood Transf. & Sn-MP Treatment');
    legend('Blood Transfusion', 'Sn-MP');
elseif(Phototherapy == 1 && SnMP == 1 && Blood_Transfusion == 0)
    title('Phototherapy. & Sn-MP Treatments');
    legend('Phototherapy', 'Sn-MP');
elseif(Blood_Transfusion == 1 && Phototherapy == 1 && SnMP == 0)
    title('Blood Transf. & Phototherapy Treatments');
    legend('Blood Transfusion', 'Phototherapy');
elseif(Blood_Transfusion == 1 && Phototherapy == 1 && SnMP == 1)
    title('Blood Transf. & Phototherapy & Sn-MP Treatments');
    legend('Blood Transfusion', 'Phototherapy', 'Sn-MP');
end

hold off
