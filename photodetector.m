clc;
close all;
clear all;

%% Data (InGaAs photodetector FGA01)
data = dlmread('InGaAs.txt');
lambda_data = data(:,1);%in nm
resp_data = data(:,2);% responsivity data in A/W

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