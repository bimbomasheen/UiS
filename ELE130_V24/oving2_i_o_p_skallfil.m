%% i)
clear; clc; close all

T =    [-150, -100, -50, -0.001, 0, 2, 4, 6, 8, 10, 30, 50, 70, 100];
rho =  [                                                           ];

plot(T,rho,'b')
xlabel('temperatur [$^{\circ}$C]')
ylabel('tetthet [kg/m$^3$]')


%% o) 

U = 0.5;
w = 0.3;

%% p)

w1 = 0.3;
w2 = 0.55;
w3 = 0.8;
f1 = w1/(2*pi);
f2 = w2/(2*pi);
f3 = w3/(2*pi);