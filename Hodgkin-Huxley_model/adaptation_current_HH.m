function [A_current_e,A_current_i] = adaptation_current_HH(spiketimes_e,spiketimes_i)
%ADAPTATION_CURRENT_HH Summary of this function goes here

% solve diff. equation for adaptation current and implement it
% parameters for adaptation current taken from:
% https://neuronaldynamics.epfl.ch/online/Ch6.S2.html


% Parameters
t_k = 10; % time constant for adaptation current, indexing over k neurons
a_k = 0; % adaptation coupling
b_k = 40; % adaptation gain
u_rest = -70; % resting potential
u_r = -55; % voltage reset




% Time series
tspan = linspace(1, 1000, 100000);
w0 = 0; % initial condition for w_k

u = @(t) 10 * sin(0.01 * t) + 1; % example: sine wave oscillating around 1


% Define and solve the differential equation:  
% t_k * dw_k / dt = a_k(u_t - u_rest) - w_k + b_k * delta(t - t_f)
odefun = @(t, w) (1/t_k) * ((a_k * (u(t) - u_rest) - w) + (t >= t_f) * b_k);
[t, w_k] = ode45(odefun, tspan, w0);

% Plot the results
figure;
plot(t, w_k, 'DisplayName', 'adaptation current');
xlabel('Time');
ylabel('Adaptation Current (w_k)');
title('Adaptation Current over Time with Varying u_t');


grid on;

% Plot u_t for reference
hold on;
plot(t, u(t), '--', 'DisplayName', 'membrane potential');
legend show;

plot(t_f, w_k(find(t >= t_f, 1)), 'ro', ...
    'MarkerSize', 8, 'DisplayName', 'spike'); % plot where delta function activates




end

