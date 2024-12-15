%% r(t) as function of I, 3 values of b_k
tic
clc
clear all

% Parameters
t_m = 20;        % Time constant for membrane potential
t_k = 100; % time constant for adaptation current, indexing over k neurons
a_k = 0.0; % adaptation coupling
u_r = -55; % voltage reset
u_rest = -70; % resting potential
u_th= -50; % spike threshold

dt=0.01; %in seconds
spk_times=[];
counter=0;
count=0;
I = zeros(1,100000);
u(1)=0; % intitial conditions u
w(1)=0;  % initial conditions adaptation current

r_start = 25000; % start of signal
r_end = 75000; % end of signal
r_width = r_end - r_start;
bin_size = 1000; % bin size in ms
max_amp = 65; % in pA, picoamperes
gain_values=(5:5:60);
rt_all = cell(1,length(gain_values));
rt_I = zeros(1,max_amp); % array of firing rates for plotting against I
R = 1.00; % Resistance, in megaohms

A = 1; % toggle for adaptation, 1 for on 0 for off


for b_k = gain_values
    amp_count = 0;
    for Amp = 0.5:0.5:max_amp
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

        firing_rate = spikes / r_width;
        rt_I(amp_count) = firing_rate;
    end
    count=count+1;
    rt_all(1,count) = {rt_I};
end


figure(1)
plot(rt_all{1,1}, 'LineWidth', 2,'LineStyle', '--');
hold on
plot(rt_all{1,2})
plot(rt_all{1,6})
plot(rt_all{1,end})
xlabel('Input Current [pA]');
ylabel('Firing rate [Hz]');
title('r(t) as a function of input current');
xlim([20 1.5*max_amp])
%ylim([-0.0005 0.004])
legend(sprintf('b_k = %g', gain_values(1)), ...
       sprintf('b_k = %g', gain_values(2)), ...
       sprintf('b_k = %g', gain_values(6)), ...
       sprintf('b_k = %g', gain_values(end)), 'Location','best'); 

toc
%% Tiled plot layout

% Create a tiled layout with 3 rows and 2 columns
tiledlayout(2, 2);

% Two smaller plots
nexttile(1);
plot(I, 'LineWidth', 1, 'LineStyle', '-'); % Adjust as needed
xlabel('time [ms]');
ylabel('amplitude [pA]');
title('Input Current');
% xlim([20 2*max_amp]);
ylim([-10 max_amp*1.2]);

nexttile(2)
plot(diff(spk_times), 'b-*')
xlim([2 158]);
ylim([0 250]);
xlabel('spike no');
title('interspike interval')

nexttile([1 2]);
plot(u, 'LineWidth', 1, 'LineStyle', '-'); % Adjust as needed
xlabel('time [ms]');
ylabel('membrane potential [mV]');
title('Voltage');
xlim([24000 42000]);
ylim([-80 -45]);
yline(-50,'--k')
legend('membrane potential','u_{th}', 'Location', 'best')


% Larger plot
% nexttile([2, 2]);
% 
% plot(rt_all{1,1}, 'LineWidth', 2, 'LineStyle', '--');
% hold on;
% plot(rt_all{1,2});
% plot(rt_all{1,6});
% plot(rt_all{1,end});
% xlabel('Input Current [pA]');
% ylabel('Firing rate [Hz]');
% title('r(t) as a function of input current');
% xlim([20 2*max_amp]);
% ylim([0 0.007]);
% legend(sprintf('b_k = %g', gain_values(1)), ...
%        sprintf('b_k = %g', gain_values(2)), ...
%        sprintf('b_k = %g', gain_values(6)), ...
%        sprintf('b_k = %g', gain_values(end)), 'Location', 'best');
