clc;
close all;
clear all;

%% Parameter

Temp = 290; %temperature in kelvin
Irr = 500; %Irradiance(Wm-2)

%% Solar cell
[Iout,Vout,Pout] = solar_func(Irr,Temp)

%% Laser
I = -Iout;
[lambda_in,Pout_laser] = laser_func(I)

%% Photodetector
[Iout_pd] = photodetector_func(lambda_in,Pout_laser,Temp)
