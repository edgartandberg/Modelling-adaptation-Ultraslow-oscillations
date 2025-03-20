% model of synapse transmission

function [I_syn_e, I_syn_i, g_e_new, g_i_new] = synaptic_current(g_e_old, g_i_old, k,  spk_times_e, spk_times_i, A)


% Parameters
g_bar_ee =  0.01; % Strength for max conductance
%g_bar_ei = 0.000001;
g_bar_ii =  0.004; % g_ee * 0.4
%g_bar_ie = 0.000001;
tau = 85 * 1e+5;         % Time constant 
E_syn_e = 0;   % pre-synaptic reversal potential for E
E_syn_i = -75; % pre-synaptic reversal potential for I
v_post = -65;  % post-synaptic membrane potential (rest)

%disp(spk_times_e)



if isempty(spk_times_e) == 1
    spk_times_e = 0;

else
    spk_times_e = spk_times_e;
end

if isempty(spk_times_i) == 1
    spk_times_i = 0;

else
    spk_times_i = spk_times_i;
end

time_points_e = spk_times_e; % Array of time points for jumps
time_points_i = spk_times_i; % Array of time points for jumps




% Initialize new conductance
g_e_new = g_e_old;
g_i_new = g_i_old;



if ~isempty(spk_times_e) && any(spk_times_e > 0)
    % Update conductance for excitatory spikes
    for t = spk_times_e
        g_e_new = g_e_new + g_bar_ee * exp(-(k-t) / tau);
    end

else 
    g_e_new = g_e_new * exp(-(k) / tau);

end

if ~isempty(spk_times_i) && any(spk_times_i > 0)
    % Update conductance for inhibitory spikes
    for t = spk_times_i
        g_i_new = g_i_new + g_bar_ii * exp(-(k-t) / tau);
    end

else
    g_i_new = g_i_new * exp(-(k) / tau);

end








% g_e = g_e(1:round(t_final/dt)); % making sure g is correct size
% g_i = g_i(1:round(t_final/dt));


I_syn_e = g_e_new .*(E_syn_e-v_post) .* A(k);
I_syn_i = g_i_new .*(E_syn_i-v_post) .* A(k);


%disp(I_syn_e)





