close all
clc
clearvars

% Parameters
t_m = 20;        % Time constant for membrane potential
theta_0 = 60; % amount threshold jumps by
t_th = 100; % threshold adaptation time constant, how sharp the decay is
u_r = -55; % voltage reset
u_rest = -70; % resting potential
dt=0.01; %in seconds
spk_times=[];
counter=0;
u(1)=40; % intitial conditions u
u_th(1)=0;  % initial conditions threshold
I(1,1:25000)=0;
I(1,25000:75000)=200;
I(1,75000:100000)=0;


uth0=50;
% Calculate v_th(t)
for t=2:100000
    u(t)= u(t-1) + dt*(-(u(t-1) - u_rest) + I(t-1))/t_m;
    %    u_th(t) = u_th(t-1) + theta_0 * exp(-(t*dt) / t_th);
    if counter == 0
        u_th(t) = uth0 + theta_0 * exp(-(t*dt) / t_th);
    else
        exp_offset=0;
        for c=1:counter
            exp_offset = exp_offset + theta_0 * exp(-(t*dt - dt*spk_times(c)) / t_th);
        end
        u_th(t) = uth0 + exp_offset;
    end
    if (u(t)>=u_th(t))
        u(t)=u_r;
        counter=counter+1;
        spk_times(counter)=t;
       
    end
end
% figure
figure();
subplot(3,1,1)
hold on;
plot(u, 'b-', 'DisplayName', 'u(t)');
plot(u_th, 'r-', 'DisplayName', 'threshold');
xlabel('Time (s)');
ylabel('');
title('Membrane potential + threshold')
legend;
grid on;
subplot(3,1,2)
plot(I, 'k-')
title('Input current');
ylim([-50 250])
subplot(3,1,3)
plot(diff(spk_times), 'b-*')
title('interspike interval')
hold off



% 
% figure;
% plot(u_th, 'r-*', 'DisplayName', 'threshold');
