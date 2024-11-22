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

inputType = 'O'; % 'O' for oscillatory input, 'S' for step input
Amplitude = 200; % amplitude of signal
r_start = 25000; % start of signal
r_end = 75000; % end of signal
r_width = r_end - r_start;

%input type
I = zeros(1,100000);
if inputType == 'S'
    I(1,r_start:r_end) = Amplitude; % step input
elseif inputType == 'O'
    t = 0:0.001:100;
    I = Amplitude * sin(2 * pi * 0.02 * t); % oscillatory input
else
    error('Invalid input type. Use ''S'' for step or ''O'' for oscillatory.');
end


% Generate noise to every nth step
noise = cos(2*pi*0.05 + 2*pi*rand) * 50 * randn(1, size(I, 2));
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
title(['Voltage w adaptation, r(t) = ', num2str(firing_rate), 'Hz'])
ylim([-90 130])
legend;
grid on;
subplot(2,2,2)
plot(I, 'k-')
title('Input current');
ylim([-50 500])
subplot(2,2,3)
plot(diff(spk_times), 'b-*')
xlabel('spike no');
title('interspike interval')
subplot(2,2,4)
plot(w, 'r-', 'DisplayName', 'w(t)');
title('adaptation')
hold off;

%% Heatmap a_k b_k

% Define ranges for a_k and b_k


dt=0.01; %in seconds
heat_spk_times=[];
heat_counter=0;
t_m = 20;        % Time constant for membrane potential
t_k = 100; % time constant for adaptation current, indexing over k neurons
u_r = -90; % voltage reset
u_rest = -70; % resting potential
a_k_values = -5.0:1:5.0; % a_k from -5.0 to 5.0
b_k_values = 0:1000:10000; % b_k from 0 to 10000
t_f = 500; % time point for spike
u_th= 50; % spike threshold

r_array = zeros(1,0);
% Loop through different values of a_k and b_k
for i = 1:length(a_k_values)
    i = a_k_values(i);
    for j = 1:length(b_k_values)
        j = b_k_values(j);
        % Initialize variables
        u(1) = 40; 
        w(1) = 0;  
        I = zeros(1, 100000); % Input current initialization
        inputType = 'S'; % 'O' for oscillatory input, 'S' for step input
        Amplitude = 200; % amplitude of signal
        r_start = 25000; % start of signal
        r_end = 75000; % end of signal
        r_width = r_end - r_start;

    %input type
    I = zeros(1,100000);
    if inputType == 'S'
        I(1,r_start:r_end) = Amplitude; % step input
    elseif inputType == 'O'
        t = 0:0.001:100;
        I = Amplitude * sin(2 * pi * 0.02 * t); % oscillatory input
    else
        error('Invalid input type. Use ''S'' for step or ''O'' for oscillatory.');
    end


        % Define and solve the differential equations:
        for t = 2:100000
            u(t) = u(t-1) + dt*(-(u(t-1) - u_rest) - w(t-1) + I(t-1))/t_m;
            if (u(t) >= u_th)
                u(t) = u_r;
                heat_counter=heat_counter+1;
                heat_spk_times(heat_counter)=t;
                
                w(t) = w(t-1) + dt*(i*(u(t-1) - u_rest) - w(t-1) + j*t_k)/t_k;
            else
                w(t) = w(t-1) + dt*(i*(u(t-1) - u_rest) - w(t-1))/t_k;
            end
            
        end
        

        %calculate firing rate
        bin_size = 1000; % bin size in ms
        spikes = length(heat_spk_times); 
        bins = r_width / bin_size;
        bin_edges = r_start:bin_size:r_end;
        heat_spike_counts = histcounts(heat_spk_times, bin_edges);
        heat_firing_rates = mean(heat_spike_counts)/ bin_size*1000; 

        % Append the current heat_firing_rates to the matrix
        r_array = [r_array heat_firing_rates];
    end

end

heatmatrix = reshape(r_array, [],length(a_k_values));

Xlabels = a_k_values;
Xlabels(mod(Xlabels,5) ~= 0) = " ";

Ylabels = b_k_values;
Ylabels(mod(Ylabels,5000) ~= 0) = " ";


h = heatmap(heatmatrix);
h.XDisplayLabels = Xlabels;
h.YDisplayLabels = Ylabels;




h.Title = 'Firing rate for varying coupling and gain';
h.XLabel = 'adaptation coupling';
h.YLabel = 'adaptation gain';
h.Colormap = autumn;

%% Heatmap t_k b_k


dt=0.01; %in seconds
heat_spk_times=[];
heat_counter=0;
t_m = 20;        % Time constant for membrane potential

u_r = -90; % voltage reset
u_rest = -70; % resting potential
a_k = 0.0; 
t_k_values = 0:10:100; 
b_k_values = 0:1000:10000; % b_k from 0 to 10000 in steps of 1000
t_f = 500; % time point for spike
u_th= 50; % spike threshold

r_array = zeros(1,0);
% Loop through different values of t_k and b_k
for i = 1:length(t_k_values)
    i = t_k_values(i);
    for j = 1:length(b_k_values)
        j = b_k_values(j);
        % Initialize variables
        u(1) = 40; 
        w(1) = 0;  
        I = zeros(1, 100000); % Input current initialization
        inputType = 'S'; % 'O' for oscillatory input, 'S' for step input
        Amplitude = 200; % amplitude of signal
        r_start = 25000; % start of signal
        r_end = 75000; % end of signal
        r_width = r_end - r_start;

    %input type
    I = zeros(1,100000);
    if inputType == 'S'
        I(1,r_start:r_end) = Amplitude; % step input
    elseif inputType == 'O'
        t = 0:0.001:100;
        I = Amplitude * sin(2 * pi * 0.02 * t); % oscillatory input
    else
        error('Invalid input type. Use ''S'' for step or ''O'' for oscillatory.');
    end


        % Define and solve the differential equations:
        for t = 2:100000
            u(t) = u(t-1) + dt*(-(u(t-1) - u_rest) - w(t-1) + I(t-1))/t_m;
            if (u(t) >= u_th)
                u(t) = u_r;
                heat_counter=heat_counter+1;
                heat_spk_times(heat_counter)=t;
                
                w(t) = w(t-1) + dt*(a_k*(u(t-1) - u_rest) - w(t-1) + j*i)/i;
            else
                w(t) = w(t-1) + dt*(a_k*(u(t-1) - u_rest) - w(t-1))/i;
            end
            
        end
        

        %calculate firing rate
        bin_size = 1000; % bin size in ms
        spikes = length(heat_spk_times); 
        bins = r_width / bin_size;
        bin_edges = r_start:bin_size:r_end;
        heat_spike_counts = histcounts(heat_spk_times, bin_edges);
        heat_firing_rates = mean(heat_spike_counts)/ bin_size*1000; 

        % Append the current heat_firing_rates to the matrix
        r_array = [r_array heat_firing_rates];
    end

end

heatmatrix = reshape(r_array, [],length(b_k_values));


h = heatmap(t_k_values,b_k_values,heatmatrix);

h.Title = 'Firing rate for varying time constant and gain';
h.XLabel = 'time constant';
h.YLabel = 'adaptation gain';
h.Colormap = autumn;