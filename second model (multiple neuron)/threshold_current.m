
% Parameters
theta_0 = -1;  % Example value for theta_0
t_th = 100;     % Example value for t_th
T = 2000;      % Total time steps

% Initialize u_th
u_th = zeros(1, T);
u_th(1) = 50; % Initial condition

% Compute u_th over time
for t = 3:T
    u_th(t) = u_th(t-1) + theta_0 * exp(-t / t_th);
    
end

% Time vector
time = 1:T;

% Plot the function
figure;
plot(time, u_th, 'b-', 'LineWidth', 2);
xlabel('Time (t)');
ylabel('u_{th}(t)');
title('Plot of u_{th}(t)');
grid on;