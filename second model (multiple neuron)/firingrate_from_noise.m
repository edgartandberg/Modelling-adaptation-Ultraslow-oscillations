%% r(t) as a function of noise

clc
clear
%tic

% Parameters
t_m = 5;        % Time constant for membrane potential
t_k = 100; % time constant for adaptation current, indexing over k neurons
a_k = 0.5; % adaptation coupling, unit is nS (nanosiemens)
b_k = 1.0; % adaptation gain
u_r = -55; % voltage reset
u_rest = -70; % resting potential
u_th= -50; % spike threshold

dt=0.01; %in seconds
spk_times=[];
counter=0;
count=0;
u(1)=0; % intitial conditions u
w(1)=0;  % initial conditions adaptation current


r_start = 20000; % start of signal
r_end = 70000; % end of signal
r_width = r_end - r_start;
bin_size = 1000; % bin size in ms
max_amp = 2.1e-6; % in pA, picoamperes
Amp_array = 0.1e-6:0.1e-6:max_amp;
gain_values=(0:5:60);
tau_values = (5:5:50);
%R = 1.00; % Resistance, in megaohms
R = 10^7; %Resistance
noise_gain_values=(0:0.1e-6:0.5e-5);
rt_all = cell(1,length(noise_gain_values));
rt_noise = zeros(1,length(Amp_array)); % array of firing rates for plotting against noise
noise_all = cell(1,length(noise_gain_values));
I_all = cell(1,length(noise_gain_values)); % all values of input current


A = 0; % toggle for adaptation, 1 for on 0 for off

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
        rt_noise(n_count) = firing_rate;
    end
    count=count+1;
    rt_all(1,count) = {rt_noise};
end


% figure(1)
plot(noise_gain_values,rt_all{1,1})
hold on
%plot(rt_all{1,2})
% plot(rt_all{1,6})
% plot(rt_all{1,13})
xlabel('Noise amplitude [A]');
ylabel('Firing rate [Hz]');
title('r(t) as a function of noise');
% xlim([0 length(noise_gain_values)])
% ylim([0 120])


legend(sprintf('b_k = %g', gain_values(1)),'Location','best'); 

% legend(sprintf('b_k = %g', gain_values(1)), ...
%        sprintf('b_k = %g', gain_values(2)), ...
%        sprintf('b_k = %g', gain_values(3)), ...
%        sprintf('b_k = %g', gain_values(4)), 'Location','best'); 



% figure(2)
% plot(I_all{1,3}, 'LineWidth', 1, 'LineStyle', '-'); % Adjust as needed
% xlabel('time [ms]');
% ylabel('amplitude [pA]');
% title('Input Current');
% % xlim([20 2*max_amp]);
% ylim([-10 max_amp*1.2]);
% 

%toc

%% Tiled layout
% 
% % Two smaller plots
% nexttile(1);
% plot(I_all{1,3}, 'LineWidth', 1, 'LineStyle', '-'); % Adjust as needed
% xlabel('time [ms]');
% ylabel('amplitude [pA]');
% title('Input Current');
% % xlim([20 2*max_amp]);
% ylim([-10 max_amp*1.2]);
% 
% nexttile(2)
% plot(diff(spk_times), 'b-*')
% xlim([2 40]);
% ylim([0 2000]);
% xlabel('spike no');
% title('interspike interval')
% 
% nexttile([1 2]);
% plot(u, 'LineWidth', 1, 'LineStyle', '-'); % Adjust as needed
% xlabel('time [ms]');
% ylabel('membrane potential [mV]');
% title('Voltage');
% xlim([20000 32000]);
% ylim([-72 -45]);
% yline(-50,'--k')
% legend('membrane potential','u_{th}', 'Location', 'best')

figure();

% Subplot for membrane potential and adaptation current
subplot(3,1,1)
hold on;
yyaxis left
plot(u, 'k-', 'DisplayName', 'Membrane potential');
xlabel('Time (ms)');
ylabel('Voltage [mV]'); % Label for left y-axis
ax = gca; % Get current axes
ax.YAxis(1).Color = 'k'; % Set right y-axis color

yyaxis right

plot(w, 'Color', '#D95319', 'DisplayName', 'threshold'); % Set color using hex code
ylabel('Threshold'); % Label for right y-axis
ax = gca; % Get current axes
ax.YAxis(2).Color = '#D95319'; % Set right y-axis color


ylabel('Adaptation current');
title('Membrane potential with adaptation current')
xlim([25000 65000]);
ylim([-80 -30])
yyaxis right
ylim([-1 12]); % Adjust limits for right y-axis if needed
legend;
grid on;

% Set limits for the right y-axis
yyaxis right
ylim([0 25]); % Adjust these values as needed for your data

% Subplot for input current
subplot(3,1,2)
plot(I, 'k-')
ylabel('Input current [pA]');
xlabel('Time (ms)');
title('Input current');
ylim([-20 82])
xlim([25000 65000]);


% Subplot for interspike interval
subplot(3,1,3)
plot(diff(spk_times), 'b-*')
title('interspike interval')
ylabel('Spike distance [ms]');
xlabel('Spike no');
ylim([200 450]);
xlim([3 60])
hold off

