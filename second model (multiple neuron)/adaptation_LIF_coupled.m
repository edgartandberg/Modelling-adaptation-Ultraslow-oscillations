% solve diff. equation for adaptation current and implement it
% parameters for adaptation current taken from:
% https://neuronaldynamics.epfl.ch/online/Ch6.S2.html

clc
clearvars

% Parameters
t_m = 20;        % Time constant for membrane potential
a_k = 0.0; % adaptation coupling
t_k = 100; % time constant for adaptation current, indexing over k neurons
b_k = 10000; % adaptation gain
u_r = -55; % voltage reset
u_rest = -70; % resting potential

t_f = 500; % time point for spike
u_th= 50; % spike threshold



dt=0.01; %in seconds
spk_times=[];
counter=0;
u(1)=u_rest; % intitial conditions u
w(1)=0;  % initial conditions adaptation current
I(1,1:25000)=0;
I(1,25000:75000)=200;
I(1,75000:100000)=0;




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

% plot
figure(2);
subplot(3,1,1)
hold on;
plot(u, 'b-', 'DisplayName', 'u(t)');
plot(w, 'r--', 'DisplayName', 'w(t)');
xlabel('Time (s)');
ylabel('');
title('Membrane potential + adaptation')
%ylim([-50 250])
legend;
grid on;
subplot(3,1,2)
plot(I, 'k-')
title('Input current');
ylim([-50 250])
subplot(3,1,3)
plot(diff(spk_times), 'b-*')
title('interspike interval')
hold off;

%% Plot for phase plane %%

% figure(2);
% plot(u, w, 'DisplayName', 'trajectory');
% xlabel('u');
% ylabel('w');
% hold on;
% plot(u(1), w(2), 'ro', 'MarkerSize', 8, 'DisplayName', 'Initial Condition'); 
% legend show;
% title('Phase Plane Plot');
% grid on;


