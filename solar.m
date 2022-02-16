function Iout = solar(Irr,Temp)
e = 1.6e-19; 
kb = 1.38e-23;

T = 300; % temperature in kelvin
I0 = 25e-9; %reverse saturation current(A)
%Iph = 10e-3; %Photocurrent of the solar cell(A)
K = 2e-5; % for Si solar cell
Iph = K*Irr; 
Eg = 1.14; % Bandgap of Si
n = 1; %Ideality factor
Rs = 10; %series resistance of the equivalent model
Rp = 1e6; %parallel resistance of the equivalent model

%% Temperature Effect
T_new = Temp;
I0_old = I0;
I0=((T_new^3)*exp(Eg./(kb*T_new/e)))*I0/((T^3)*exp(Eg/(kb*T/e))); %reverse saturation current(A)

%% Calculation of current, power(considering Rp and Rs)
V=0:0.001:0.4;
I_total = zeros(1,length(V));
for i = 1:length(V)
    fcn = @(I) -I - Iph + I0*(exp(e*(V(i)-I*Rs)/(n*kb*T))-1) + (V(i)-I*Rs)/Rp;
    I = fzero(fcn,Iph);
    I_total(i)= I;
end
Power = (-I_total.*V);

%% Load Line
R = 30;
err = (-V/R-I_total);
index = find(abs(err)<0.05e-3);

%% I-V Curve Plot

plot(V,I_total*1e3,'Linewidth',2)
xlabel('Voltage, V(V)')
ylabel('Current,I_{total}(mA)')
grid on;
hold on
line([V(1), V(end)], [0, 0], 'Color', [0,0,0],'LineStyle','-.','linewidth',2);
plot(V,(-V/R)*1e3);
plot(V(index),I_total(index)*1e3,'ro')


Iout  = I_total(index)*1e3; % in mA

end