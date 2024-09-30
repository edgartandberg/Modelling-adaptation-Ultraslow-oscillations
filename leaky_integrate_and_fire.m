
%tic

clear all
close all
clc

% Play with the input current, add noise

tau_m = 5; % time constant in ms
t_ref = 10; % refractory period in ms, change for adjusting spike frequency
N = 10000; % number of measurements
Iext = zeros(1,N); % initialize membrane current
I_0 = 1.6e-6; % input current
R = 10^7; % resistance in ohms
a_gain = t_ref*2; % adaptation gain, set to *0 for no adaptation


u_rest = -70 ;% resting membrane potential in mV
u_th = -55 ;% spiking threshold
u_hp = -90 ;% hyperpolarization
u_spike = 20 ;% membrane potential at spike peak

t = linspace(-0.2*N, 0.8*N, N); % time series




%%

noise = 5*1e-8*randn(size(Iext)); % noise when current is injected

 

num_start = 0.2 * (N); % start of input current
num_end   = 0.8 * (N);  % end of input current
num_width = round(num_end - num_start);

num_pulses = 3; % number of pulses
length_pulse = round(num_width / num_pulses); % only in use for step input

freq = 0.001 * num_pulses; % adjusting frequency for right amount of pulses

time = num_start:1:num_end;
signal =  I_0*sin(time*freq);


Iext(num_start : num_start + num_width) = signal; % modelling input current as sine function for oscillations
% change to = I_0 for step input

%%

Iext = Iext + noise;
Iext(1:num_start) = 0; % resetting signal before and after current to 0

% add or remove 50 - 53 for without/with oscillatory input

% Iext(num_start + length_pulse : num_start + 2 * length_pulse) = 0;
% Iext(num_start + 3 * length_pulse : num_start + 4 * length_pulse) = 0;
% Iext(num_start + 5 * length_pulse : num_start + 6 * length_pulse) = 0;
% Iext(num_start + 7 * length_pulse : num_start + 8 * length_pulse) = 0;
Iext(num_end:end) = 0;


RI = R.*Iext;




%% LiF model as ODE

%du/dt = (-u + u_rest + R*Iext)/tau_m

u = u_rest + RI.*(1-exp(-(t/tau_m))); %solution for above equation for u at u(0) = u_rest

% Loops over time steps, spikes if threshold is reached, before
% hyperpolarizing and returning to resting potential


a = t_ref;
for ii = 1 : num_end - 2 - a % Adjust loop range to avoid out-of-bounds errors
    if u(ii + 1) > u_th
        u(ii) = u_spike;
        u(ii + 1) = u_hp;
        u(ii + 2 : ii + 2 + a) = u_rest;
        a = round(a + a_gain); % Increment by adaptation gain
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

%toc