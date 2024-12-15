%% r(t) as a function of noise

clc
clear
%tic

% Parameters
t_m = 20;        % Time constant for membrane potential
t_k = 100; % time constant for adaptation current, indexing over k neurons
a_k = 0.5; % adaptation coupling, unit is nS (nanosiemens)
b_k = 5.0; % adaptation gain
u_r = -55; % voltage reset
u_rest = -70; % resting potential
u_th= -50; % spike threshold

dt=0.01; %in seconds
spk_times=[];
counter=0;
count=0;
u(1)=0; % intitial conditions u
w(1)=0;  % initial conditions adaptation current

r_start = 25000; % start of signal
r_end = 75000; % end of signal
r_width = r_end - r_start;
bin_size = 1000; % bin size in ms
max_amp = 65;
gain_values=(0:5:60);
noise_gain_values=(0:0.1:30);
rt_all = cell(1,length(noise_gain_values));
rt_noise = zeros(1,max_amp); % array of firing rates for plotting against noise
R = 1.0; % Resistance, in megaohms
noise_all = cell(1,length(noise_gain_values));
I_all = cell(1,length(noise_gain_values)); % all values of input current

% For the example of adaptation with noise: 
% amp = 37, noise_gain_values=(0:0.1:0.4)

% For example of firing rate over noise amp, for different adaptation
% gains:
% amp = 65, noise_gain_values=(0:0.1:30)


%input type
I = zeros(1,100000);
n = 0;
I(1,r_start:r_end) = max_amp;

for b_k = gain_values
    n_count = 0;
    for n = noise_gain_values
        
        noise = n * randn(1, size(I, 2));
        
        for idx = 1:1:length(I)
            I(idx) = I(idx) + noise(idx);
        end
        I_all(1,n_count+1) = {I}; 
        %disp(I)
        noise_all(1,n_count+1) = {noise};
        n_count = n_count +1;
        clear spk_times
        counter = 0;

        for t=2:100000
            u(t)= u(t-1) + dt*(-(u(t-1) - u_rest) - w(t-1) + R.*I(t-1))/t_m;
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
        % Firing rate   
        spikes = length(spk_times); % number of spikes
        firing_rate = spikes / r_width;
        rt_noise(n_count) = firing_rate;
    end
    count=count+1;
    rt_all(1,count) = {rt_noise};
end


% figure(1)
plot(rt_all{1,1})
hold on
%plot(rt_all{1,2})
plot(rt_all{1,6})
plot(rt_all{1,13})
xlabel('Noise amplitude');
ylabel('Firing rate [Hz]');
title('r(t) as a function of noise');
% xlim([0 length(noise_gain_values)])
% ylim([0 120])
legend(sprintf('b_k = %g', gain_values(1)), ...
       sprintf('b_k = %g', gain_values(2)), ...
       sprintf('b_k = %g', gain_values(3)), ...
       sprintf('b_k = %g', gain_values(4)), 'Location','best'); 



figure(2)
plot(I_all{1,3}, 'LineWidth', 1, 'LineStyle', '-'); % Adjust as needed
xlabel('time [ms]');
ylabel('amplitude [pA]');
title('Input Current');
% xlim([20 2*max_amp]);
ylim([-10 max_amp*1.2]);


%toc

%% Tiled layout

% Two smaller plots
nexttile(1);
plot(I_all{1,3}, 'LineWidth', 1, 'LineStyle', '-'); % Adjust as needed
xlabel('time [ms]');
ylabel('amplitude [pA]');
title('Input Current');
% xlim([20 2*max_amp]);
ylim([-10 max_amp*1.2]);

nexttile(2)
plot(diff(spk_times), 'b-*')
xlim([2 40]);
ylim([0 2000]);
xlabel('spike no');
title('interspike interval')

nexttile([1 2]);
plot(u, 'LineWidth', 1, 'LineStyle', '-'); % Adjust as needed
xlabel('time [ms]');
ylabel('membrane potential [mV]');
title('Voltage');
xlim([20000 32000]);
ylim([-72 -45]);
yline(-50,'--k')
legend('membrane potential','u_{th}', 'Location', 'best')