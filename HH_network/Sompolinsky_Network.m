%% 3rd pass at network. Changed some cell arrays to arrays
clc
clearvars

tic


% 1000ms sim takes ~1 hour
t_final=550; % in ms
dt=0.01;
dt05=dt/2;
m_steps=round(t_final/dt);


max_amps = 15.0; % 
 
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
% 
% visualizeConnectivityMatrix(connectivity)


%%


v_e = zeros(1, m_steps);
v_i = zeros(1, m_steps);
v_e(1)=-70; % initial condition
v_i(1)=-70; % initial condition

% spiketrains_e_all = cell(1, size_network); % initialize cell arrays for neurons
% spiketrains_i_all = cell(1, size_network);

spiketrains_e_all = zeros(size_network,m_steps);
spiketrains_i_all = zeros(size_network,m_steps);



for i = 1:size_network

    spiketrains_e_all(i,:) = v_e;
    spiketrains_i_all(i,:) = v_i;

end

% % initialize arrays for gating variables, excitatory
m_e = zeros(1,m_steps);
h_e = zeros(1,m_steps);
n_e = zeros(1,m_steps);
z_e = zeros(1,m_steps);
b_e = zeros(1,m_steps);

% inhibitory
m_i = zeros(1,m_steps);
h_i = zeros(1,m_steps);
n_i = zeros(1,m_steps);
z_i = zeros(1,m_steps);
b_i = zeros(1,m_steps);

% cell arrays to save gating variables for each neuron
m_e_all = cell(1, size_network);
h_e_all = cell(1, size_network);
n_e_all = cell(1, size_network);
z_e_all = cell(1, size_network);
b_e_all = cell(1, size_network);

m_i_all = cell(1, size_network);
h_i_all = cell(1, size_network);
n_i_all = cell(1, size_network);
z_i_all = cell(1, size_network);
b_i_all = cell(1, size_network);


for i = 1:size_network

    m_e_all{i} = m_e;
    h_e_all{i} = h_e;
    n_e_all{i} = n_e;
    z_e_all{i} = z_e;
    b_e_all{i} = b_e;
    
    m_i_all{i} = m_i;
    h_i_all{i} = h_i;
    n_i_all{i} = n_i;
    z_i_all{i} = z_i;
    b_i_all{i} = b_i;


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

A = ones(1, m_steps); %Normalization constant



%%
clc
tic

i_ext_e = max_amps;

z(1)=0.0;
a(1)=0.0;
b(1)=0.0;
s(1)=0.0;

z_e = z;
z_i = z;

b_e = b;
b_i = b;


for k=2:m_steps
    for  j = 1:size_network
    i_ext_e = max_amps;
    %i_ext_i = max_amps;

    v_e = spiketrains_e_all(j,:);
    v_i = spiketrains_i_all(j,:);

    m_e = m_e_all{j};
    h_e = h_e_all{j};
    n_e = n_e_all{j};
    z_e = z_e_all{j};
    b_e = b_e_all{j};

    m_i = m_i_all{j};
    h_i = h_i_all{j};
    n_i = n_i_all{j};
    z_i = z_i_all{j};
    b_i = b_i_all{j};

    

    
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
    
    
    % I_syn_e = 0;
    % for i = 1:length(connections_synaptic)
    %     I_syn_e = I_syn_e + I_syn_e_all(connections_synaptic(i),k);
    %     disp(i)
    % end


    % how long external current lasts, in ms

current_length = 5000;


    if k > current_length
        i_ext_e = 0;
        i_ext_i = 0;
    end



    if j > 1 % only first neuron gets ext input
        i_syn_e = I_syn_e + size_network*I_syn_i;
    else
        i_syn_e = i_ext_e + I_syn_e + size_network*I_syn_i;
    end

i_syn_i = i_ext_i + 3*I_syn_e + size_network*I_syn_i; 
% 
% disp(i_syn_e)
% disp(i_syn_i)


% Without adaptation
    % [v_e, t_e, e_counter, m_e, h_e, n_e, spk_times_e] = excitatory_HH_pass2(v_e, k, t_final,dt,i_syn_e, m_e, h_e, n_e);
    % [v_i, t_i, i_counter, m_i, h_i, n_i, spk_times_i] = inhibitory_HH_pass2(v_i, k, t_final,dt,i_syn_i, m_i, h_i, n_i);




% With adaptation
    [v_e, t_e, e_counter, m_e, h_e, n_e, spk_times_e, i_ext_e_all,I_z, z_e, b_e] = E_Sompolinsky_adaptation(v_e, k,t_final,dt,i_syn_e, m_e, h_e, n_e, z_e, b_e);
    [v_i, t_i, i_counter, m_i, h_i, n_i, spk_times_i, i_ext_i_all, z_i, b_i] = I_Sompolinsky_adaptation(v_i, k, t_final,dt,i_syn_i, m_i, h_i, n_i, z_i, b_i);

    %disp(spk_times_e)


    %generate synaptic currents

    g_e = g_e_all{j};
    g_i = g_i_all{j};

    
    [I_syn_e, I_syn_i, g_e_new, g_i_new] = synaptic_current(v_e, v_i, g_e, g_i, k, spk_times_e, spk_times_i, A); % generating synaptic current

    %disp(I_syn_e)


    g_e = g_e_new;
    g_i = g_i_new;


    g_e_all{j} = g_e;
    g_i_all{j} = g_i;


    I_syn_e_all(j, k) = I_syn_e; % Store the synaptic current for excitatory
    I_syn_i_all(j, k) = I_syn_i; % Store the inhibitory synaptic current
    


    spiketrains_e_all(j,1:length(v_e)) = v_e;
    spiketrains_i_all(j,1:length(v_i)) = v_i;

    i_ext_e_all(j,k) = i_ext_e;
    i_ext_i_all(j,k) = i_ext_i;


   % gating variable, commented out for performance
    m_e_all{j} = m_e;
    h_e_all{j} = h_e;
    n_e_all{j} = n_e;
    b_e_all{j} = b_e;
    z_e_all{j} = z_e;

    m_i_all{j} = m_i;
    h_i_all{j} = h_i;
    n_i_all{j} = n_i;
    b_i_all{j} = b_i;
    z_i_all{j} = z_i;

    end


end

t_e = t_e(1:end-1);


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



%% Calculate sequence duration for excitatory raster
% Checks where the first sequence starts in the first neuron, then takes 
% time difference from first spike in first neuron to last spike in last neuron.
% Spikes are binned into 'bursts' for each first and last neuron, to
% calculate the sequence duration.
% 8000 in line 325 here is 80 ms in the simulation.


spikes_bursts_e = [];

for i = 1:length(spike_times_e{1, 1})-1
    if (spike_times_e{1, 1}(i+1)) - (spike_times_e{1, 1}(i)) <= 4000 
        spikes_bursts_e = [spikes_bursts_e, i];
    end
end

% Initialize variables
bursts_e = {}; % Cell array to hold bursts
current_burst_e = []; % Array to hold current burst

% Loop through the array
for i = 1:length(spikes_bursts_e)
    % Add current number to the current burst
    current_burst_e(end + 1) = spikes_bursts_e(i);

    % Check if the next number is more than 1 larger than the current
    if i < length(spikes_bursts_e) && (spikes_bursts_e(i + 1) - spikes_bursts_e(i) > 1)
        % Append the last value + 1 to the current burst
        current_burst_e(end + 1) = current_burst_e(end) + 1;

        % Store the current burst and reset
        bursts_e{end + 1} = current_burst_e; 
        current_burst_e = []; % Reset for the next burst
    end
end

if ~isempty(current_burst_e)
    current_burst_e(end + 1) = current_burst_e(end) + 1; % Append last value + 1
    bursts_e{end + 1} = current_burst_e;
end


% --------------------------------------------- %
% Binning spikes in last neuron in bursts as well
% --------------------------------------------- %



last_bursts_e = [];
    
% Bursts analysis for last neuron as well
for i = 1:length(spike_times_e{size_network, 1})-1
    if (spike_times_e{size_network, 1}(i+1)) - (spike_times_e{size_network, 1}(i)) <= 4000  
        last_bursts_e = [last_bursts_e, i];
    end
end

% Initialize variables
final_bursts_e = {}; 
current_burst_e = []; 

% Loop through the array
for i = 1:length(last_bursts_e)
    % Add current number to the current burst
    current_burst_e(end + 1) = last_bursts_e(i);

    % Check if the next number is more than 1 larger than the current
    if i < length(last_bursts_e) && (last_bursts_e(i + 1) - last_bursts_e(i) > 1)
        % Append the last value + 1 to the current burst
        current_burst_e(end + 1) = current_burst_e(end) + 1;

        % Store the current burst and reset
        final_bursts_e{end + 1} = current_burst_e; 
        current_burst_e = []; % Reset for the next burst
    end
end

if ~isempty(current_burst_e)
    current_burst_e(end + 1) = current_burst_e(end) + 1; % Append last value + 1
    final_bursts_e{end + 1} = current_burst_e;
end


time_window = 8000; % time window to check next spike in sequence

% Array to store the indices of neuron 1 spikes that start a sequence
sequence_starts_e = [];
last_spike_indices_e = [];

% --------------------------------------------- %
% Calculating durations from binned bursts
% --------------------------------------------- %

results_sequences_e = [];

for i = 1:length(bursts_e)
    % Check if neuron 1 has spikes in the current burst
    if ~isempty(bursts_e{i})
        valid_sequence_count = 0; % Counter for valid sequences
        first_spike_time = []; % To store the first spike time
        last_spike_time = []; % To store the last spike time

        % Loop through each spike in neuron 1 in all bursts
        for spike_idx = 1:length(bursts_e{i})
            spike_time = spike_times_e{1}(bursts_e{i}(spike_idx));
            valid_sequence = true; % Flag to check if the sequence is valid
            last_spike_index_e = [];

            % Check spikes in subsequent neurons
            for neuron_idx = 2:length(spike_times_e)
                next_spikes = spike_times_e{neuron_idx};
                % Check if there are spikes in the next neuron within the time window
                if isempty(next_spikes) || all(next_spikes < spike_time | next_spikes > spike_time + time_window)
                    valid_sequence = false; % No valid spike found in this neuron
                    % disp('neuron')
                    % disp(neuron_idx)
                    % disp('burst')
                    % disp(i)
                    % disp('spike')
                    % disp(spike_idx)
                    % disp('spike')
                    % disp(spike_time)
                    break; % Exit the loop if the sequence is broken
                end
                % Update the spike_time for the next neuron
                spike_time = next_spikes(find(next_spikes >= spike_time, 1)); % Get the first valid spike
                if neuron_idx == length(spike_times_e)
                    last_spike_index_e = find(next_spikes == spike_time, 1); % Get the index of the spike in the last neuron
                end
            end

            % If a valid sequence was found, increment the valid sequence count
            if valid_sequence
                valid_sequence_count = valid_sequence_count + 1; % Increment count
                if spike_idx == 1
                    first_spike_time = spike_times_e{1}(bursts_e{i}(1)); % Save first spike time
                end
                if neuron_idx == length(spike_times_e)
                    % last_spike_time = spike_time; % Save last spike time
            
                    last_spike_time = spike_times_e{size_network}(final_bursts_e{i}(end));
                end
            end
        end

        % Check if the number of valid sequences is equal to or greater than the number of spikes in the burst
        if valid_sequence_count >= length(bursts_e{i})

            last_spike_time = round(last_spike_time*dt);
            first_spike_time = round(first_spike_time*dt);
            time_difference = last_spike_time - first_spike_time; % Calculate time difference
            % Save the results (burst index, first spike time, last spike time, time difference)
            results_sequences_e(end + 1, :) = [i, first_spike_time, last_spike_time, time_difference]; % Store results
        end
    end

end

%% duration for inhibitory raster

spikes_bursts_i = [];

for i = 1:length(spike_times_i{1, 1})-1
    if (spike_times_i{1, 1}(i+1)) - (spike_times_i{1, 1}(i)) <= 4000 
        spikes_bursts_i = [spikes_bursts_i, i];
    end
end

% Initialize variables
bursts_i = {}; % Cell array to hold bursts
current_burst_i = []; % Array to hold current burst

% Loop through the array
for i = 1:length(spikes_bursts_i)
    % Add current number to the current burst
    current_burst_i(end + 1) = spikes_bursts_i(i);

    % Check if the next number is more than 1 larger than the current
    if i < length(spikes_bursts_i) && (spikes_bursts_i(i + 1) - spikes_bursts_i(i) > 1)
        % Append the last value + 1 to the current burst
        current_burst_i(end + 1) = current_burst_i(end) + 1;

        % Store the current burst and reset
        bursts_i{end + 1} = current_burst_i; 
        current_burst_i = []; % Reset for the next burst
    end
end

if ~isempty(current_burst_i)
    current_burst_i(end + 1) = current_burst_i(end) + 1; % Append last value + 1
    bursts_i{end + 1} = current_burst_i;
end


% --------------------------------------------- %
% Binning spikes in last neuron in bursts as well
% --------------------------------------------- %



last_bursts_i = [];
    
% Bursts analysis for last neuron as well
for i = 1:length(spike_times_i{size_network, 1})-1
    if (spike_times_i{size_network, 1}(i+1)) - (spike_times_i{size_network, 1}(i)) <= 4000  
        last_bursts_i = [last_bursts_i, i];
    end
end

% Initialize variables
final_bursts_i = {}; 
current_burst_i = []; 

% Loop through the array
for i = 1:length(last_bursts_i)
    % Add current number to the current burst
    current_burst_i(end + 1) = last_bursts_i(i);

    % Check if the next number is more than 1 larger than the current
    if i < length(last_bursts_i) && (last_bursts_i(i + 1) - last_bursts_i(i) > 1)
        % Append the last value + 1 to the current burst
        current_burst_i(end + 1) = current_burst_i(end) + 1;

        % Store the current burst and reset
        final_bursts_i{end + 1} = current_burst_i; 
        current_burst_i = []; % Reset for the next burst
    end
end

if ~isempty(current_burst_i)
    current_burst_i(end + 1) = current_burst_i(end) + 1; % Append last value + 1
    final_bursts_i{end + 1} = current_burst_i;
end


time_window = 8000; % time window to check next spike in sequence

% Array to store the indices of neuron 1 spikes that start a sequence
sequence_starts_i = [];
last_spike_indices_i = [];




results_sequences_i = [];

for i = 1:length(bursts_i)
    % Check if neuron 1 has spikes in the current burst
    if ~isempty(bursts_i{i})
        valid_sequence_count = 0; % Counter for valid sequences
        first_spike_time = []; % To store the first spike time
        last_spike_time = []; % To store the last spike time

        % Loop through each spike in neuron 1 in all bursts
        for spike_idx = 1:length(bursts_i{i})
            spike_time = spike_times_i{1}(bursts_i{i}(spike_idx));
            valid_sequence = true; % Flag to check if the sequence is valid
            last_spike_index_i = [];

            % Check spikes in subsequent neurons
            for neuron_idx = 2:length(spike_times_e)
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

            % If a valid sequence was found, increment the valid sequence count
            if valid_sequence
                valid_sequence_count = valid_sequence_count + 1; % Increment count
                if spike_idx == 1
                    first_spike_time = spike_times_i{1}(bursts_i{i}(1)); % Save first spike time
                end
                if neuron_idx == length(spike_times_i)
                    % last_spike_time = spike_time; % Save last spike time
            
                    last_spike_time = spike_times_i{size_network}(final_bursts_i{i}(end));
                end
            end
        end

        % Check if the number of valid sequences is equal to or greater than the number of spikes in the burst
        if valid_sequence_count >= length(bursts_i{i})

            last_spike_time = round(last_spike_time*dt);
            first_spike_time = round(first_spike_time*dt);
            time_difference = last_spike_time - first_spike_time; % Calculate time difference
            % Save the results (burst index, first spike time, last spike time, time difference)
            results_sequences_i(end + 1, :) = [i, first_spike_time, last_spike_time, time_difference]; % Store results
        end
    end

end


%% Plot for excitatory raster

figure()

plot(results_sequences_i(:,1), results_sequences_i(:,end), 'k', 'LineWidth', 3); % Plot in black with line width of 3
xlabel('Sequence#'); % Label for x-axis
ylabel('Sequence Length [ms]'); % Label for y-axis
title('Sequence Durations as Function of Sequence Index'); % Title of the plot
ylim([0, results_sequences_i(end, end) * 1.25]); 


% Set x-axis ticks to whole number values from results(:,1)
xticks(unique(round(results_sequences_i(:,1)))); % Use unique whole number values

% Adjust x-axis limits to add space
xlim([min(results_sequences_i(:,1)) - 1, max(results_sequences_i(:,1)) + 1]); % Add space to the left and right

toc