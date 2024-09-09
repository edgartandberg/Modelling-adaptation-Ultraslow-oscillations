
clear
close all
clc


tau_m = 10; % time constant in ms
t_ref = 5; % refractory period in ms
N = 500; % number of measurements
I_0 = 2e-6; % input current
Iext = zeros(1,N); % initialize membrane current
R = 10^7 ;% resistance in ohms


u_rest = -70 ;% resting membrane potential in mV
u_th = -40 ;% spiking threshold
u_hp = -90 ;% hyperpolarization
u_spike = 20 ;% membrane potential at spike peak

t = linspace(-200, 800, N); % time series




%%

num_start = 0.2 * N;
num_end   = 0.22 * N;  
num_width = round(num_end - num_start);
Iext(num_start : num_start + num_width) = I_0;

RI = R.*Iext;



%% LiF model as ODE

%du/dt = (-u + u_rest + R*Iext)/tau_m

u = u_rest + RI.*(1-exp(-(t/tau_m))); %solution for above equation for u at u(0) = u_rest


%% Plotting input current and voltage

figure(1)
subplot(2,1,1);
plot(t,Iext,'r','LineWidth',3)

title('Input current')
xlabel('time')
ylabel('I_0')
ylim([0 I_0])
yline(0)

subplot(2,1,2);
plot(t,u, 'b', 'LineWidth',2)
title('Voltage')
xlabel('time')
ylabel('u(t)')
ylim([u_hp u_spike])