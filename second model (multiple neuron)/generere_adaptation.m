%tic

clearvars
close all
clc

% Play with the input current, add noise

tau_m = 5; % time constant in ms
t_ref = 7; % refractory period in ms, change for adjusting spike frequency
N = 10000; % number of measurements
Iext = zeros(1,N); % initialize membrane current
I_0 = 1.6e-6; % input current
R = 10^7; % resistance in ohms
a_gain = t_ref*3; % adaptation gain, set to *0 for no adaptation


u_rest = -70 ;% resting membrane potential in mV
u_th = -55 ;% spiking threshold
u_hp = -90 ;% hyperpolarization
u_spike = 20 ;% membrane potential at spike peak

t = linspace(-0.2*N, 0.8*N, N); % time series

spiketrains = 25; % number of spiketrains we want to simulate. Change to 1 for single neuron
num_start = 0.2 * (N); % start of input current
num_end   = 0.8 * (N);  % end of input current
num_pulses = 5; % number of pulses

%% Building connectivity matrix

figure(1)
plotConnectivityMatrix(spiketrains);

%% Generate input current + Function for leaky integrate and fire + adaptation, LIF.m

% Initialize output variables
RI_all = cell(spiketrains, 1);
Iext_all = cell(spiketrains, 1);
u_all = cell(spiketrains, 1);
a_all = cell(spiketrains, 1);
spikecount_all = cell(spiketrains, 1);

for i = 1:spiketrains
    [RI, Iext] = Gen_Current(R, Iext, num_start, num_end, num_pulses, I_0);
    
    RI_all{i} = RI;
    Iext_all{i} = Iext;

    [u, a, spikecount] = LIF(u_rest, RI, t, tau_m, t_ref, num_end, u_th, u_spike, u_hp, a_gain, num_start);

    u_all{i} = u;
    a_all{i} = a;
    spikecount_all{i} = spikecount;
end

%% Plotting input current and voltage
figure(2) 
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
hold off


%% Raster plot

figure(3)
xlim([-0.1*N num_end*1.05]);
ylim([0 spiketrains])

timeStepS = 1;
rasterPlot(spikecount_all, timeStepS);

%toc