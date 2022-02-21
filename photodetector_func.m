function [Iout] = photodetector_func(lambda_in,Pout_laser,Temp)
%% Data (InGaAs photodetector FGA01)
data = dlmread('InGaAs.txt');
lambda_data = data(:,1);%in nm
resp_data = data(:,2);% responsivity data in A/W
figure
plot(lambda_data,resp_data,'Linewidth',2)
xlabel('Lambda(nm)')
ylabel('Responsivity(A/W)')
title('Spectral Response of Photodetector')
%% parameter
e = 1.6e-19; 
kb = 1.38e-23;
n = 1; %Ideality factor
T = 300; % temperature in kelvin
I0 = 25e-9; %reverse saturation current(A)
Eg = 0.784; %In0.5Ga0.5As
%% Temperature Effect
T_new = Temp;
I0_old = I0;
I0=((T_new^3)*exp(Eg./(kb*T_new/e)))*I0_old/((T^3)*exp(Eg/(kb*T/e))); %reverse saturation current(A)


%% Iph calculation
lambda_in = lambda_in*1e9; % from laser (in nm)
%Pout_laser = 0.0012;% from laser (in W)

err = lambda_in - lambda_data;
index = find(err == min(abs(err)));
R = resp_data(index);

Iph = R*Pout_laser;

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
figure
plot(V,I_total*1e3,'Linewidth',2)
xlabel('Voltage, V(V)')
ylabel('Current,I_{total}(mA)')
grid on;
hold on
line([V(1), V(end)], [0, 0], 'Color', [0,0,0],'LineStyle','-.','linewidth',2);
plot(V(index),I_total(index)*1e3,'ro')
title('I-V characteristics of Photodetector')
Iout = I_total(index); % in A

end