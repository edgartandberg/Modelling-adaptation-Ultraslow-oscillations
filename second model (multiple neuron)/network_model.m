%tic

clearvars
close all
clc


t_m = 5; % time constant in ms
t_ref = 50; % refractory period in ms, change for adjusting spike frequency
N = 800; % number of measurements
Iext = zeros(1,N); % initialize membrane current
I_0 = 2.0e-6; % input current
R = 10^7; % resistance in ohms
%R = 1.0;
a_gain = 0*t_ref*5; % set to *0 for no adaptation and *1 for adaptation


u_rest = -70 ;% resting membrane potential in mV
u_th = -50 ;% spiking threshold
u_hp = -90 ;% hyperpolarization
u_spike = 20 ;% membrane potential at spike peak



t = linspace(-0.2*N, 0.8*N, N); % time series

spiketrains = 1; % number of spiketrains we want to simulate. Change to 1 for single neuron
num_start = 200; % start of input current
num_end   = 700;  % end of input current
num_pulses = 2; % number of pulses

%% Building connectivity matrix

% figure(1)
% [connectivitymatrix] = plotConnectivityMatrix(spiketrains);


%% Generate input current + Function for leaky integrate and fire + adaptation, LIF.m

% Initialize output variables
RI_all = cell(spiketrains, 1);
Iext_all = cell(spiketrains, 1);
u_all = cell(spiketrains, 1);
a_all = cell(spiketrains, 1);
spikecount_all = cell(spiketrains, 1);

for i = 1:spiketrains
    [RI, Iext, noise] = Gen_Current(R, Iext, num_start, num_end, num_pulses, I_0);

    RI_all{i} = RI;
    Iext_all{i} = Iext;

    [u, a, spikecount] = LIF(u_rest, R, Iext, t_m, t_ref, num_end, u_th, u_spike, u_hp);
    u_all{i} = u;
    a_all{i} = a;
    spikecount_all{i} = spikecount;
end

%% Switch to this part if you want network + input current with delay on each spike train
% at 0 for no delay

delay_increment = 50; % Adjust this value as needed

for i = 1:spiketrains

    % Calculate the new start time with delay
    current_start = num_start + (i-1)*delay_increment; 

    % current for 1st spike train
    [RI, Iext, noise] = Gen_Current(R, Iext, num_start, num_end, num_pulses, I_0);

    RI_all{i} = RI;
    Iext_all{i} = Iext;

    if i > 1
        if connectivitymatrix(i,i-1) > 0.90

            previous_chain = u_all{i-1,1};

            weight = connectivitymatrix(i,i-1); 

            [u,a,spikecount] = Synaptic_input(previous_chain,u,t_ref,u_rest,num_end,u_th,u_spike,current_start,u_hp,a_gain,weight,num_start);
        else
            u = 0;

        end
    else

    % membrane potential for 1st spiketrain
    [u, a, spikecount] = LIF(u_rest, R, Iext, t_m, t_ref, num_end, u_th, u_spike, u_hp);


    end

    u_all{i} = u;
    a_all{i} = a;
    spikecount_all{i} = spikecount;
end

%% Plotting input current and voltage

figure(2) 
subplot(2,1,1);
plot(t,Iext_all{1,1},'k','LineWidth',3)
%title('Input current')
xlabel('Time [ms]')
ylabel('Input current [A]')
xlim([-0.1*N num_end*1.05])
ylim([0 2.5e-6])
yline(0)

subplot(2,1,2);
plot(t,u_all{1,1}, 'k', 'LineWidth',2)
%title('Voltage')
xlabel('Time [ms]')
ylabel('V [mV]')
xlim([-0.1*N num_end*1.05])
ylim([u_hp-10 u_spike+5])
yline(u_th,'--k','Threshold')
hold off


%% Raster plot

% figure(3)
% xlim([-0.1*N num_end*1.05]);
% ylim([0 spiketrains])
% 
% timeStepS = 1;
% rasterPlot(spikecount_all, timeStepS);

%toc


