
%% 2nd pass at network
clc
clearvars

tic

t_final=200; % in ms
dt=0.01;
dt05=dt/2;
m_steps=round(t_final/dt);


max_amps = 7; % in

i_ext_e = max_amps*ones(1,m_steps);
i_ext_i = 0;



size_network = 6; % amount of both excitatory and inhibitory neurons, size of network
i_ext_e_all = zeros(size_network,m_steps);
i_ext_i_all = zeros(size_network,m_steps);

%%

% Make connectivity
n = size_network; % Number of neurons
connectivity = createConnectionMatrix(n);
%disp(connectivity);

%visualizeConnectivityMatrix(connectivity)


%%


v_e = zeros(1, m_steps);
v_i = zeros(1, m_steps);
v_e(1)=-70; % initial condition
v_i(1)=-70; % initial condition

spiketrains_e_all = cell(1, size_network); % initialize cell arrays for neurons
spiketrains_i_all = cell(1, size_network);

for i = 1:size_network

    spiketrains_e_all{i} = v_e;
    spiketrains_i_all{i} = v_i;

end

% % initialize arrays for gating variables, excitatory
m_e = zeros(1,m_steps);
h_e = zeros(1,m_steps);
n_e = zeros(1,m_steps);

% inhibitory
m_i = zeros(1,m_steps);
h_i = zeros(1,m_steps);
n_i = zeros(1,m_steps);

% cell arrays to save gating variables for each neuron
m_e_all = cell(1, size_network);
h_e_all = cell(1, size_network);
n_e_all = cell(1, size_network);

m_i_all = cell(1, size_network);
h_i_all = cell(1, size_network);
n_i_all = cell(1, size_network);

for i = 1:size_network

    m_e_all{i} = m_e;
    h_e_all{i} = h_e;
    n_e_all{i} = n_e;
    
    
    m_i_all{i} = m_i;
    h_i_all{i} = h_i;
    n_i_all{i} = n_i;

end



% intitialize arrays for conductances, excitatory
g_e_all = cell(1, size_network);
I_syn_e_all = zeros(size_network, m_steps);

% inhibitory
g_i_all = cell(1, size_network);
I_syn_i_all = zeros(size_network, m_steps);


% initital conditions for gating variables, excitatory
m_e(1)=alpha_m(v_e(1))/(alpha_m(v_e(1))+beta_m(v_e(1)));
h_e(1)=0.5; 
n_e(1)=0.5;

% inhibitory
m_i(1)=alpha_m(v_i(1))/(alpha_m(v_i(1))+beta_m(v_i(1)));
h_i(1)=0.5; 
n_i(1)=0.5; 

g_e = zeros(1, m_steps); % intialize synaptic conductance
g_i = zeros(1, m_steps);


for i = 1:size_network

    g_e_all{i} = g_e;
    g_i_all{i} = g_i;

end


I_syn_e = zeros(1, m_steps);
I_syn_i = zeros(1, m_steps);

A = 10*ones(1, m_steps); %Normalization constant



%%

for i = 1:size_network
    if i == size_network
        i_ext_e = max_amps;

    else
        i_ext_e = 0;

    end
end



for k=2:m_steps
    for  j = 1:size_network
    % i_ext_e = max_amps;
    % i_ext_i = max_amps;

    v_e = spiketrains_e_all{j};

    m_e = m_e_all{j};
    h_e = h_e_all{j};
    n_e = n_e_all{j};

    m_i = m_i_all{j};
    h_i = h_i_all{j};
    n_i = n_i_all{j};




    
    % [v_e, t_e, e_counter, m_e, h_e, n_e, spk_times_e] = excitatory_HH_pass2(v_e, k, t_final,dt,i_ext_e, m_e, h_e, n_e);
    % [v_i, t_i, i_counter, m_i, h_i, n_i, spk_times_i] = inhibitory_HH_pass2(v_i, k, t_final,dt,i_ext_i, m_i, h_i, n_i);
    



    % run them all through excitatory_HH and inhibitory_HH again with the
    % external current + synaptic current, for every neuron.

    connections_synaptic = [];

    for i = 1:size_network
        if connectivity(j,i) == 1   
           connections_synaptic = [connections_synaptic i];
        elseif connectivity(j,i) == 0
           connections_synaptic = connections_synaptic;
        end
    end
    
    %disp(connections_synaptic)
    
    % for i = 1:length(connections_synaptic)
    %     I_syn_e = I_syn_e + I_syn_e_all(connections_synaptic(i),k-1);
    % end
    
    I_syn_e = 0;
    for i = 1:length(connections_synaptic)
        I_syn_e = I_syn_e + I_syn_e_all(connections_synaptic(i),k-1);
    end


    % how long external current lasts, in ms

current_length = 500;


    if k > current_length
        i_ext_e = 0;
    end

    % i_syn_e = i_ext_e + I_syn_e + size_network*I_syn_i;

    if j > 1 % only first neuron gets ext input
        i_syn_e = I_syn_e + size_network*I_syn_i;
    else
        i_syn_e = i_ext_e + I_syn_e + size_network*I_syn_i;
    end

    % if j == 1
    %     i_syn_e = i_ext_e + I_syn_e + size_network*I_syn_i;
    % else
    %     i_syn_e = I_syn_e + size_network*I_syn_i;
    % 
    % end


    %i_syn_e = i_ext_e + I_syn_e + size_network * I_syn_i;  % correct one from before

    %i_syn_e = i_ext_e + I_syn_e;

    
    
    [v_e, t_e, e_counter, m_e, h_e, n_e, spk_times_e] = excitatory_HH_pass2(v_e, k, t_final,dt,i_syn_e, m_e, h_e, n_e);


    % inhibitory neurons

    i_syn_i = i_ext_i + 3*I_syn_e + size_network*I_syn_i; 
    [v_i, t_i, i_counter, m_i, h_i, n_i, spk_times_i] = inhibitory_HH_pass2(v_i, k, t_final,dt,i_syn_i, m_i, h_i, n_i);

      % generate synaptic currents


    g_e_old = g_e_all{j};
    g_i_old = g_i_all{j};

    g_e_old = g_e_old(k-1);
    g_i_old = g_i_old(k-1);
    
    [I_syn_e, I_syn_i, g_e_new, g_i_new] = synaptic_current(v_e, v_i, g_e_old, g_i_old, k, spk_times_e, spk_times_i, A); % generating synaptic current
    
  
    g_e(k) = g_e_new;
    g_i(k) = g_i_new;


    g_e_all{j} = g_e;
    g_i_all{j} = g_i;
    
    
    I_syn_e_all(j, k) = I_syn_e; % Store the synaptic current for excitatory
    I_syn_i_all(j, k) = I_syn_i; % Store the inhibitory synaptic current


    % I_syn_e_all(j, k) = I_syn_e; % Store the synaptic current for each network
    % I_syn_i_all(j, k) = I_syn_i; % Store the inhibitory synaptic current 
    spiketrains_e_all{j} = v_e;
    spiketrains_i_all{j} = v_i;

    i_ext_e_all(j,k) = i_ext_e;
    i_ext_i_all(j,k) = i_ext_i;


   % gating variable
    m_e_all{j} = m_e;
    h_e_all{j} = h_e;
    n_e_all{j} = n_e;

    m_i_all{j} = m_i;
    h_i_all{j} = h_i;
    n_i_all{j} = n_i;

    end


end


%% Plot

[spike_times_e, spike_times_i] = rasterPlot_HH(spiketrains_e_all,spiketrains_i_all, current_length ...
    );


% figure(3)
% plot(t_e, v_e, 'LineWidth', 2);
% xlabel('Time (ms)');
% ylabel('Amplitude');
% title('post-synaptic conductance g, from pre-synaptic spike train');
% grid on;
% %axis([0 0.1 0 45]); % Axis limits

toc

%%

% %pause(10)
% 
% % Create a figure
% f = figure;
% 
% % Set the figure size (width, height)
% f.Position = [500, 200, 800, 700]; % [left, bottom, width, height]
% % Loop through each subplot
% for i = 1:6
%     subplot(6, 1, 7 - i); % Create a subplot in a 6x1 grid
%     plot(spiketrains_e_all{1, i}); % Plot the data
%     title(['Spike Train ', num2str(i)]); % Add a title for each subplot
%     xlabel('Time'); % Label for x-axis
%     ylabel('Amplitude'); % Label for y-axis
% end
% 
% % Adjust layout for better visibility
% %tight_layout(); % Optional: Use this if you have the 'tight_layout' function available
%% Time between spikes and velocity

% Time window for checking spikes
time_window = 2000;

% Array to store the indices of neuron 1 spikes that start a sequence
sequence_starts_e = [];
last_spike_indices_e = [];


% Check if neuron 1 has spikes
if ~isempty(spike_times_e{1})
    % Loop through each spike in neuron 1
    for spike_idx = 1:length(spike_times_e{1})
        spike_time = spike_times_e{1}(spike_idx);
        valid_sequence = true; % Flag to check if the sequence is valid

        last_spike_index_e = [];

        % Check spikes in subsequent neurons
        for neuron_idx = 2:length(spike_times_e)
            next_spikes = spike_times_e{neuron_idx};
            % Check if there are spikes in the next neuron within the time window
            if isempty(next_spikes) || all(next_spikes < spike_time | next_spikes > spike_time + time_window)
                valid_sequence = false; % No valid spike found in this neuron
                break; % Exit the loop if the sequence is broken
            end
            % Update the spike_time for the next neuron
            spike_time = next_spikes(find(next_spikes >= spike_time, 1)); % Get the first valid spike
            if neuron_idx == length(spike_times_e)
                last_spike_index_e = find(next_spikes == spike_time, 1); % Get the index of the spike in the last neuron
            end
        end

        % If a valid sequence was found, save the index of the spike in neuron 1
        if valid_sequence
            sequence_starts_e(end + 1) = spike_idx; % Append the index
            last_spike_indices_e(end + 1) = last_spike_index_e; % Append the index of the last neuron
        end
    end
end

% Sequence duration for I neurons

sequence_starts_i = [];
last_spike_indices_i = [];


% Check if neuron 1 has spikes
if ~isempty(spike_times_i{1})
    % Loop through each spike in neuron 1
    for spike_idx = 1:length(spike_times_i{1})
        spike_time = spike_times_i{1}(spike_idx);
        valid_sequence = true; % Flag to check if the sequence is valid

        last_spike_index_i = [];

        % Check spikes in subsequent neurons
        for neuron_idx = 2:length(spike_times_i)
            next_spikes = spike_times_i{neuron_idx};
            % Check if there are spikes in the next neuron within the time window
            if isempty(next_spikes) || all(next_spikes < spike_time | next_spikes > spike_time + time_window)
                valid_sequence = false; % No valid spike found in this neuron
                break; % Exit the loop if the sequence is broken
            end
            % Update the spike_time for the next neuron
            spike_time = next_spikes(find(next_spikes >= spike_time, 1)); % Get the first valid spike
            if neuron_idx == length(spike_times_i)
                last_spike_index_i = find(next_spikes == spike_time, 1); % Get the index of the spike in the last neuron
            end
        end

        % If a valid sequence was found, save the index of the spike in neuron 1
        if valid_sequence
            sequence_starts_i(end + 1) = spike_idx; % Append the index
            last_spike_indices_i(end + 1) = last_spike_index_i; % Append the index of the last neuron
        end
    end
end


% disp('Indices of neuron 1 spikes that start a sequence:');
% disp(sequence_starts);


%% Sequence durations

% 
%  % choose which sequences based on spikes in neuron 1
% for i = 1:length(sequence_starts_e)
%     spike_diff = spike_times_e{1, size_network}(last_spike_indices_e(i)) - spike_times_e{1, 1}(sequence_starts_e(i));
%     spike_diffs_e(1,sequence_starts_e(i)) = spike_diff;
% end
% 
% 
% 
% M_all  = mean(spike_diffs_e);
% 
% 
% figure()
% pointSize = 100; % Size of the points
% pointColor = 'b'; % Color of the points (can also use RGB triplet or hex code)
% 
% % Create the scatter plot with specified size and color
% scatter(sequence_starts_e, spike_diffs_e, pointSize, pointColor, 'filled');
% title('Sequence duration, for each sequence in simulation'); % Add a title for each subplot
% %text(1,M_true*1.04, ['mean duration = ' num2str(M_true)])
% yline(M_all, 'LineWidth', 2, 'LineStyle','--','Color','k')
% xlabel('Sequence duration by sequence index'); % Label for x-axis
% ylabel('Time Difference [ms]'); % Label for y-axis
% 
% velocity = (spike_times_e{1, size_network}(end) - spike_times_e{1,1}(end-1));

