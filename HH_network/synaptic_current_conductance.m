% model of synapse transmission

function [I_syn_e, I_syn_i, g_e_new, g_i_new] = synaptic_current_conductance(g_bar_ee, g_bar_ii, g_e, g_i, k,  spk_times_e, spk_times_i, A)


% Parameters
% g_bar_ee =  0.2; % Strength for max conductance
% %g_bar_ei = 0.000001;
% g_bar_ii =  0.08; % g_ee * 0.4
% %g_bar_ie = 0.000001;
tau = 1000;         % Time constant, *0.01 for tau in ms
E_syn_e = 0;   % pre-synaptic reversal potential for E
E_syn_i = -75; % pre-synaptic reversal potential for I
v_post = -65;  % post-synaptic membrane potential (rest)


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



% Initialize new conductance
g_e_new = g_e;
g_i_new = g_i;


if ~isempty(spk_times_e) && any(spk_times_e)
    % Update conductance for excitatory spikes
    for t = spk_times_e(spk_times_e <= k)
        g_e_new(k) = g_e_new(k) + g_bar_ee * exp(-(k - t) / tau);
    end
elseif k > 1
    % Exponential decay when no spikes
    g_e_new(k) = g_e_new(k - 1) * exp(-1 / tau);
end

if ~isempty(spk_times_i) && any(spk_times_i)
    % Update conductance for inhibitory spikes
    for t = spk_times_i(spk_times_i <= k)
        g_i_new(k) = g_i_new(k) + g_bar_ii * exp(-(k - t) / tau);
    end
elseif k > 1
    % Exponential decay when no spikes
    g_i_new(k) = g_i_new(k - 1) * exp(-1 / tau);
end




I_syn_e = g_e_new(k) .*(E_syn_e-v_post) .* A(k);
I_syn_i = g_i_new(k) .*(E_syn_i-v_post) .* A(k);

