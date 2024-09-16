
clear all
close all
clc

% Play with the input current, add noise

tau_m = 5; % time constant in ms
t_ref = 10; % refractory period in ms, change for adjusting spike frequency
N = 5000; % number of measurements
Iext = zeros(1,N); % initialize membrane current
I_0 = 1.4e-6; % input current
R = 10^7; % resistance in ohms
a_gain = 1.3; % adaptation gain, set to 1 for no adaptation


u_rest = -70 ;% resting membrane potential in mV
u_th = -55 ;% spiking threshold
u_hp = -90 ;% hyperpolarization
u_spike = 20 ;% membrane potential at spike peak

t = linspace(-0.2*N, 0.8*N, N); % time series




%%

noise = 5*1e-8*randn(size(Iext)); % noise when current is injected

num_start = 0.20 * N; % start of input current
num_end   = 0.6* N;  % end of input current
num_width = round(num_end - num_start);
Iext(num_start : num_start + num_width) = I_0;
Iext = Iext + noise;
Iext(1:num_start) = 0; % resetting signal before and after current to 0
Iext(num_end:end) = 0;


RI = R.*Iext;




%% LiF model as ODE

%du/dt = (-u + u_rest + R*Iext)/tau_m

u = u_rest + RI.*(1-exp(-(t/tau_m))); %solution for above equation for u at u(0) = u_rest

% Loops over time steps, spikes if threshold is reached, before
% hyperpolarizing and returning to resting potential


a = t_ref;
for ii = 1 : num_end
    if u(ii+1) > u_th
       u(ii) = u_spike;
       u(ii+1) = u_hp;
       u(ii+2:ii+2+a) = u_rest;
       a = round(a*a_gain); % incrementing by adaptation gain
    end
end



%% Plotting input current and voltage

figure(1)
subplot(2,1,1);
plot(t,Iext,'r','LineWidth',3)


title('Input current')
xlabel('time')
ylabel('I_ext')
xlim([-0.1*N num_end*1.05])
ylim([0 I_0*1.25])
yline(0)


subplot(2,1,2);
plot(t,u, 'b', 'LineWidth',2)
title('Voltage')
xlabel('time')
ylabel('u(t)')
xlim([-0.1*N num_end*1.05])
ylim([u_hp-10 u_spike+5])
yline(u_th,'--k','Threshold')

