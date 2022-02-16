clc;
close all;
clear all;

%% Data (InGaAs photodetector FGA01)
data = dlmread('InGaAs.txt');
lambda_data = data(:,1);%in nm
resp_data = data(:,2);% responsivity data in A/W

%% parameter
e = 1.6e-19; 
kb = 1.38e-23;
n = 1; %Ideality factor
T = 300; % temperature in kelvin
I0 = 25e-9; %reverse saturation current(A)
%% Iph calculation
lambda_in = (1.57e-6)*1e9; % from laser (in nm)
Pout = 0.0012;% from laser (in W)

err = lambda_in - lambda_data;
index = find(err == min(abs(err)));
R = resp_data(index);

Iph = R*Pout;

%% Output voltage
RL = 30;

Vout = Iph*RL;

%% Calculation of current, power
Vr = 1.5;
V = -3:0.01:0.35;
I_total = -Iph + I0.*(exp(e*V/(n*kb*T))-1);
Power = (-I_total.*V);
index = find(V == -Vr);
%% I-V Curve Plot

plot(V,I_total*1e3,'Linewidth',2)
xlabel('Voltage, V(V)')
ylabel('Current,I_{total}(mA)')
grid on;
hold on
line([V(1), V(end)], [0, 0], 'Color', [0,0,0],'LineStyle','-.','linewidth',2);
plot(V(index),I_total(index)*1e3,'ro')