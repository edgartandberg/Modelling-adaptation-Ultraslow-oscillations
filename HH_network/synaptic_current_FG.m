% model of synapse transmission

function [I_syn_e, g_e] = synaptic_current_FG(v_e, v_i,t_final,dt, spk_times_e, spk_times_i, g_bar)
v_e = v_e(1:end); % adjust size of voltage arrays if needed
v_i = v_i(1:end);

% disp(spk_times_e)
% disp(spk_times_i)

% Parameters

        % Strength for max conductance
tau = 50;         % Time constant 
t = 0:dt:t_final;     % Time vector from 0 to t_final seconds with dt ms steps
E_syn_e = 0;
E_syn_i = -75;
v_post = -65;

time_points_e = t(spk_times_e); % Array of time points for jumps
time_points_i = t(spk_times_i); % Array of time points for jumps

% Create Heaviside step functions for each time point
H_e = sum(heaviside(t - time_points_e'), 1); % Sum of Heaviside functions
%H_i = sum(heaviside(t - time_points_i'), 1); % Same for inhibitory



% Initialize conductance arrays
g_e = zeros(size(t));
%g_i = zeros(size(t));


% Calculate conductance with decay only after jumps, excitatory
for i = 1:length(time_points_e)
    jump_time = time_points_e(i);   
    start_index = find(t >= jump_time, 1);

    % Apply decay from the activation of heaviside function
    g_e(start_index:end) = g_e(start_index:end) + g_bar * exp(-(t(start_index:end) - jump_time) / tau);
end

% Inhibitory
% for i = 1:length(time_points_i)
%     jump_time = time_points_i(i);
%     start_index = find(t >= jump_time, 1);
% 
%     % Apply decay from the activation of heaviside function
%     g_i(start_index:end) = g_i(start_index:end) + g_bar * exp(-(t(start_index:end) - jump_time) / tau);
% end


g_e = g_e(1:round(t_final/dt)); % making sure g is correct size
%g_i = g_i(1:round(t_final/dt));



I_syn_e = g_e .*(E_syn_e-v_post);
%I_syn_i = g_i .*(E_syn_i-v_post);





