% solve diff. equation for adaptation current and implement it
% parameters for adaptation current taken from:
% https://neuronaldynamics.epfl.ch/online/Ch6.S2.html

clc
clearvars

% Parameters
t_m = 20;        % Time constant for membrane potential
a_k = 0.0; % adaptation coupling
t_k = 100; % time constant for adaptation current, indexing over k neurons
b_k = 5000; % adaptation gain
u_r = -90; % voltage reset
u_rest = -70; % resting potential

t_f = 500; % time point for spike
u_th= 50; % spike threshold



dt=0.01; %in seconds
spk_times=[];
counter=0;
u(1)=40; % intitial conditions u
w(1)=0;  % initial conditions adaptation current

inputType = 'S'; % Change to 'O' for oscillatory input, 'S' for step input
Amplitude = 150; % amplitude of signal
r_start = 25000; % start of signal
r_end = 75000; % end of signal
r_width = r_end - r_start;

%input type
I = zeros(1,100000);
if inputType == 'S'
    I(1,r_start:r_end) = Amplitude; % step input
elseif inputType == 'O'
    t = 0:0.001:100;
    I = Amplitude * sin(2 * pi * 0.05 * t); % oscillatory input
else
    error('Invalid input type. Use ''S'' for step or ''O'' for oscillatory.');
end


% Generate noise to every nth step
noise = cos(2*pi*0.05 + 2*pi*rand) * 10 * randn(1, size(I, 2));
n = 1; % interval between every noisy input 

for idx = 1:n:length(I)
    I(idx) = I(idx) + noise(idx);
end


% Define and solve the differential equations:

for t=2:100000
    u(t)= u(t-1) + dt*(-(u(t-1) - u_rest) - w(t-1) + I(t-1))/t_m;
    if (u(t)>=u_th)
        u(t)=u_r;
        counter=counter+1;
        spk_times(counter)=t;
        w(t)= w(t-1) + dt*(a_k*(u(t-1) - u_rest) -w(t-1) + b_k*t_k)/t_k;
    else
        %w(t)= w(t-1) + dt*(a_k*(u(t-1) - u_rest) -w(t-1) + b_k*t_k*length(spk_times))/t_k;
        w(t)= w(t-1) + dt*(a_k*(u(t-1) - u_rest) -w(t-1))/t_k;
    end
end

%% Firing rate

bin_size = 1000; % bin size in ms
spikes = length(spk_times); % number of spikes
bins = r_width / bin_size;
bin_edges = r_start:bin_size:r_end;
spike_counts = histcounts(spk_times, bin_edges);
firing_rate = mean(spike_counts)/ bin_size*1000; 


%% Plot
figure();
subplot(2,2,1)
hold on;
plot(u, 'b-', 'DisplayName', 'u(t)');
xlabel('Time (ms)');
ylabel('');
title(['Voltage w adaptation, firing rate = ', num2str(firing_rate), 'Hz'])
ylim([-90 130])
legend;
grid on;
subplot(2,2,2)
plot(I, 'k-')
title('Input current');
ylim([-50 250])
subplot(2,2,3)
plot(diff(spk_times), 'b-*')
xlabel('spike no');
title('interspike interval')
subplot(2,2,4)
plot(w, 'r-', 'DisplayName', 'w(t)');
title('adaptation')
hold off;

