%% r(t) as function of I, 3 values of b_k
tic
clc
clear all

% Parameters
t_k = 100; % time constant for adaptation current, indexing over k neurons
a_k = 0.0; %5e-12; % adaptation coupling, e-12 = pA
u_r = -55; % voltage reset
u_rest = -70; % resting potential
u_th= -50; % spike threshold
u_spike = 20; % voltage at spike
u_hp = -90 ;% hyperpolarization



dt=0.01; %in seconds
spk_times=[];
counter=0;
count=0;
I = zeros(1,100000);
u(1)=0; % intitial conditions u
w(1)=0;  % initial conditions adaptation current

r_start = 20000; % start of signal
r_end = 70000; % end of signal
r_width = r_end - r_start;
bin_size = 1000; % bin size in ms
max_amp = 3.0e-6; % in muA, e-6 = micro amps
Amp_array = 0:1e-6:max_amp;
gain_values=(0:5e-6:1e-5);
tau_values = (5:5:50);
rt_all = cell(1,length(gain_values));
rt_I = zeros(1,length(Amp_array)); % array of firing rates for plotting against I
%R = 1.00; % Resistance, in megaohms
R = 10^7; %Resistance

A = 1; % toggle for adaptation, 1 for on 0 for off

%% Looping over gain values

t_m = 5;        % Time constant for membrane potential


for b_k = gain_values
    amp_count = 0;
    for Amp = Amp_array
        amp_count = amp_count +1 ;
        clear spk_times
        I(1,r_start:r_end) = Amp;
        counter = 0;
        for t=2:100000
            u(t) = u(t-1) + dt*(-(u(t-1) - u_rest) - A.*R.*w(t-1) + R.*I(t-1))/t_m;
            
            if (u(t)>=u_th)
                u(t)=u_r;
                counter=counter+1;
                spk_times(counter)=t;
                w(t)= A*(w(t-1) + dt*(a_k*(u(t-1) - u_rest) -w(t-1) + b_k*t_k)/t_k);
           else
              %w(t)= w(t-1) + dt*(a_k*(u(t-1) - u_rest) -w(t-1) + b_k*t_k*length(spk_times))/t_k;
              w(t)= A*(w(t-1) + dt*(a_k*(u(t-1) - u_rest) -w(t-1))/t_k);
           end
        end

        % Firing rate   
        spikes = length(spk_times); % number of spikes

        firing_rate = spikes / (r_width*dt*0.01);
        rt_I(amp_count) = firing_rate;
    end
    count=count+1;
    rt_all(1,count) = {rt_I};
end


%% Looping over tau values

b_k = 15;

for t_m = tau_values
    amp_count = 0;
    for Amp = 0.5:0.5:max_amp
        amp_count = amp_count +1 ;
        clear spk_times
        I(1,r_start:r_end) = Amp;
        counter = 0;
        for t=2:100000
            u(t) = u(t-1) + dt*(-(u(t-1) - u_rest) - A.*R.*w(t-1) + R.*I(t-1))/t_m;
            
            if (u(t)>=u_th)
                u(t) = u_r;
                counter=counter+1;
                spk_times(counter)=t;
                w(t)= A*(w(t-1) + dt*(a_k*(u(t-1) - u_rest) -w(t-1) + b_k*t_k)/t_k);
           else
              %w(t)= w(t-1) + dt*(a_k*(u(t-1) - u_rest) -w(t-1) + b_k*t_k*length(spk_times))/t_k;
              w(t)= A*(w(t-1) + dt*(a_k*(u(t-1) - u_rest) -w(t-1))/t_k);
           end
        end

        % Firing rate   
        spikes = length(spk_times); % number of spikes

        firing_rate = spikes / r_width;
        rt_I(amp_count) = firing_rate;
    end
    count=count+1;
    rt_all(1,count) = {rt_I};
end




%% FI curve plot for adaptation values

figure(1)
plot(Amp_array,rt_all{1,1}, 'LineWidth', 2, 'LineStyle', '-', 'Marker', '*', 'MarkerSize', 3);
hold on
plot(Amp_array,rt_all{1,2}, 'LineWidth', 2, 'LineStyle', '-', 'Marker', '*', 'MarkerSize', 3);
plot(Amp_array,rt_all{1,6}, 'LineWidth', 2, 'LineStyle', '-', 'Marker', '*', 'MarkerSize', 3);
plot(Amp_array,rt_all{1,51}, 'LineWidth', 2, 'LineStyle', '-', 'Marker', '*', 'MarkerSize', 3);

xlabel('Input Current [A]');
ylabel('Firing rate [Hz]');
title('Firing rate as a function of input current');
%xlim([20 1.5*max_amp])
%ylim([-0.0005 0.004])

%legend(sprintf('b_k = %g', gain_values(1)), 'Location','best'); 


legend(sprintf('b_k = %g', gain_values(1)), ...
       sprintf('b_k = %g', gain_values(2)), ...
       sprintf('b_k = %g', gain_values(6)), ...
       sprintf('b_k = %g', gain_values(51)), 'Location','best'); 

toc

%% FI curve plot for tau values


figure(1)
plot(rt_all{1,1}, 'LineWidth', 2,'LineStyle', '--');
hold on
plot(rt_all{1,2})
plot(rt_all{1,5})
plot(rt_all{1,end})
xlabel('Input Current [pA]');
ylabel('Firing rate [Hz]');
title('r(t) as a function of input current');
xlim([20 1.5*max_amp])
%ylim([-0.0005 0.004])
legend(sprintf('t_m = %g', tau_values(5)), ...
       sprintf('t_m = %g', tau_values(12)), ...
       sprintf('t_m = %g', tau_values(20)), ...
       sprintf('t_m = %g', tau_values(end)), 'Location','best'); 

%% Tiled plot layout

figure();

% Subplot for membrane potential and adaptation current
subplot(3,1,1)
hold on;
yyaxis left
plot(u, 'k-', 'DisplayName', 'Membrane potential');
xlabel('Time (ms)');
ylabel('Voltage [mV]'); % Label for left y-axis
ylim([-80 -48])
ax = gca; % Get current axes
ax.YAxis(1).Color = 'k'; % Set right y-axis color

yyaxis right
plot(w, 'Color', '#D95319', 'DisplayName', 'Adaptation current, b_k = 1 \muA'); % Set color using hex code
ylabel('Adaptation current'); % Label for right y-axis
ax = gca; % Get current axes
ax.YAxis(2).Color = '#D95319'; % Set right y-axis color


ylabel('Adaptation current');
title('Membrane potential with adaptation current')
xlim([15000 75000]);
yyaxis right
ylim([-10e-8 5e-8]); % Adjust limits for right y-axis if needed
legend;
grid on;

% Subplot for input current
subplot(3,1,2)
plot(I, 'k-')
ylabel('Input current [A]');
xlabel('Time (ms)');
title('Input current');
ylim([-2e-6 3e-6])
xlim([15000 75000]);


% Subplot for interspike interval
subplot(3,1,3)
plot(diff(spk_times), 'b-*')
title('interspike interval')
ylabel('Spike distance [ms]');
xlabel('Spike no');
ylim([640 850]);
xlim([3 60])
hold off

