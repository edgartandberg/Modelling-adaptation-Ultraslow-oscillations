
%% 2nd pass at network - sequence duration over multiple amps
clc
clearvars

tic

t_final = 300; % in ms
dt=0.01;
dt05=dt/2;
m_steps=round(t_final/dt);

max_amps = 20; % in
i_ext_e = max_amps*ones(1,m_steps);

start_conductance = 0.04;
end_conductance   = 0.12;
step_size = 0.01;


i_ext_i = 0;

size_network = 6; % amount of both excitatory and inhibitory neurons, size of network
i_ext_e_all = zeros(size_network,m_steps);
i_ext_i_all = zeros(size_network,m_steps);

%%

% Make connectivity
n = size_network; % Number of neurons
connectivity = createConnectionMatrix(n);


%visualizeConnectivityMatrix(connectivity)


%%


v_e = zeros(1, m_steps);
v_i = zeros(1, m_steps);
v_e(1)=-70; % initial condition
v_i(1)=-70; % initial condition

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

A = 1*ones(1, m_steps); %Normalization constant



%%

conductance_array = (start_conductance:step_size:end_conductance);

spike_times_e_all = cell(length(conductance_array),size_network);
spike_times_i_all = cell(length(conductance_array),size_network);



% Looping over array of conductance
for ii = 1:length(conductance_array)
    g_bar_ii = conductance_array(ii);
    %g_bar_ii = conductance_array(ii);
    g_bar_ee = 0.2;

    %disp(i_ext_e)
    i_ext_e = max_amps;

    for k=2:m_steps
        for  j = 1:size_network
        
        v_e = spiketrains_e_all(j,:);
        v_i = spiketrains_i_all(j,:);
    
        m_e = m_e_all{j};
        h_e = h_e_all{j};
        n_e = n_e_all{j};
    
        m_i = m_i_all{j};
        h_i = h_i_all{j};
        n_i = n_i_all{j};
    
    
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

        % if v_e ~= v_i
        %     disp('true')
        % end
        
        [v_e, t_e, e_counter, m_e, h_e, n_e, spk_times_e] = excitatory_HH_pass2(v_e, k, t_final,dt,i_syn_e, m_e, h_e, n_e);
    
    
        % inhibitory neurons
        i_syn_i = i_ext_i + 3*I_syn_e + size_network*I_syn_i; 
        [v_i, t_i, i_counter, m_i, h_i, n_i, spk_times_i] = inhibitory_HH_pass2(v_i, k, t_final,dt,i_syn_i, m_i, h_i, n_i);

        % generate synaptic currents
    
    
        g_e = g_e_all{j};
        g_i = g_i_all{j};
    


        [I_syn_e, I_syn_i, g_e_new, g_i_new] = synaptic_current_conductance(g_bar_ee, g_bar_ii, g_e, g_i, k,  spk_times_e, spk_times_i, A); % generating synaptic current
        



      
        g_e = g_e_new;
        g_i = g_i_new;

    
    
        g_e_all{j} = g_e;
        g_i_all{j} = g_i;
        
        
        I_syn_e_all(j, k) = I_syn_e; % Store the synaptic current for excitatory
        I_syn_i_all(j, k) = I_syn_i; % Store the inhibitory synaptic current
    
    
        % I_syn_e_all(j, k) = I_syn_e; % Store the synaptic current for each network
        % I_syn_i_all(j, k) = I_syn_i; % Store the inhibitory synaptic current 
        spiketrains_e_all(j,:) = v_e;
        spiketrains_i_all(j,:) = v_i;


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
        for j = 1:size_network
        % Process excitatory spiketrains
        spike_times_e_all{ii,j} = detect_spikes(spiketrains_e_all(j,:));

        % Process inhibitory spiketrains
        spike_times_i_all{ii,j} = detect_spikes(spiketrains_i_all(j,:));
        end
end


%% Time between spikes and velocity. Need to calculate for all amps.

% Time window for checking spikes
time_window = 2000;

% initiate cell arrays for sequence index in first
% and last neuron. Start in column 1, end in 2.
% amps along rows

seq_e_points_I = cell(length(conductance_array),2);
seq_i_points_I = cell(length(conductance_array),2);

for a = 1:length(conductance_array)

% Array to store the indices of neuron 1 spikes that start a sequence
sequence_starts_e = [];
last_spike_indices_e = [];


% Check if neuron 1 has spikes
if ~isempty(spike_times_e_all{a})
    % Loop through each spike in neuron 1
    for spike_idx = 1:length(spike_times_e_all{a})
        spike_time = spike_times_e_all{a}(spike_idx);
        valid_sequence = true; % Flag to check if the sequence is valid

        last_spike_index_e = [];

        % Check spikes in subsequent neurons
        for neuron_idx = 2:size_network
            next_spikes = spike_times_e_all{a,neuron_idx};
            % Check if there are spikes in the next neuron within the time window
            if isempty(next_spikes) || all(next_spikes < spike_time | next_spikes > spike_time + time_window)
                valid_sequence = false; % No valid spike found in this neuron
                break; % Exit the loop if the sequence is broken
            end
            % Update the spike_time for the next neuron
            spike_time = next_spikes(find(next_spikes >= spike_time, 1)); % Get the first valid spike
            if neuron_idx == size_network
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
if ~isempty(spike_times_i_all{a})
    % Loop through each spike in neuron 1
    for spike_idx = 1:length(spike_times_i_all{a})
        spike_time = spike_times_i_all{a}(spike_idx);
        valid_sequence = true; % Flag to check if the sequence is valid

        last_spike_index_i = [];

        % Check spikes in subsequent neurons
        for neuron_idx = 2:size_network
            next_spikes = spike_times_i_all{a,neuron_idx};
            % Check if there are spikes in the next neuron within the time window
            if isempty(next_spikes) || all(next_spikes < spike_time | next_spikes > spike_time + time_window)
                valid_sequence = false; % No valid spike found in this neuron
                break; % Exit the loop if the sequence is broken
            end
            % Update the spike_time for the next neuron
            spike_time = next_spikes(find(next_spikes >= spike_time, 1)); % Get the first valid spike
            if neuron_idx == size_network
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



seq_e_points_I{a,1} = sequence_starts_e;
seq_e_points_I{a,2} = last_spike_indices_e;

seq_i_points_I{a,1} = sequence_starts_i;
seq_i_points_I{a,2} = last_spike_indices_i;


end


%% Calculate duration

for ii = 1:length(conductance_array)
    for i = 1:length(seq_e_points_I{ii})
        spike_diff_e = spike_times_e_all{ii, size_network}(seq_e_points_I{ii,2}(i)) - spike_times_e_all{ii, 1}(seq_e_points_I{ii,1}(i));
        spike_diffs_e_all{ii,seq_e_points_I{ii}(i)} = spike_diff_e;
    end
end

for ii = 1:length(conductance_array)
    for i = 1:length(seq_i_points_I{ii})
        spike_diff_i = spike_times_i_all{ii, size_network}(seq_i_points_I{ii,2}(i)) - spike_times_i_all{ii, 1}(seq_i_points_I{ii,1}(i));
        spike_diffs_i_all{ii,seq_i_points_I{ii}(i)} = spike_diff_i;
    end
end

spike_diffs_i_all = cellfun(@(x) x * 0.01, spike_diffs_i_all, 'UniformOutput', false);


%% Scatter plot for durations, I neurons



%Plot inhibitory neurons
figure;
hold on;

%color = [0 1 1];
[numRows, numCols] = size(spike_diffs_i_all);

% Initialize a vector to hold the average values
averages_i = zeros(numRows, 1);

% Loop through each row
for row = 1:numRows
    % Extract the values from the current row
    yValues = cell2mat(spike_diffs_i_all(row, :));

    % Calculate the average for the current row
    averages_i(row) = mean(yValues);

    % Initialize colors for the current row
    colors = zeros(length(yValues), 3); % m-by-3 matrix for current row

    % Calculate colors for each point in the current row
    for col = 1:length(yValues)
        % Create a gradient from magenta to blue
        colors(col, :) = [(1 - col / (length(yValues) - 1)), 0, (col / (length(yValues) - 1))]; % RGB triplet for magenta to blue
    end

    % Create a scatter plot for the current row
    scatter(conductance_array(row) * ones(size(yValues)), yValues, 100, colors, 'filled');
    hold on; % Keep the current plot
end

colors(end,:) = colors(end-1,:);

    average_array = [];

for i = 1:length(spike_diffs_i_all(:,1))
    average_array = [average_array conductance_array(i)];

end


% Plot the average line
%avgLine_i = plot(average_array, averages_i, 'k-', 'LineWidth', 2, 'DisplayName', 'Mean sequence duration');


% Change font properties
ax = gca; % Get current axes
ax.FontName = 'Helvetica-Narrow'; % Set font name
ax.FontSize = 14; % Set font size

% Create a color bar
clim([1 numCols]); % Set color axis limits
colormap(colors); % Set the colormap to the calculated colors
cb = colorbar; % Create the color bar

% Set the color bar labels to indicate the first and last columns
cb.Ticks = [1, numCols]; % Set ticks at the start and end
cb.TickLabels = {'Earlier sequence', sprintf('Later sequence', numCols)}; % Label the ticks


hold off;
grid on;
xlabel('Conductance [mS/cm^2]');
ylabel('Sequence duration [ms]');
xlim([conductance_array(1)-step_size conductance_array(end)+step_size]); 
title('Sequence durations as function of excitatory conductance, inhibitory population');
xticks(conductance_array(1):step_size:conductance_array(end));

%legend(avgLine_i); % This will display the legend with 'Averages'


%% Scatter plot for durations, E neurons


%Plot inhibitory neurons
figure;
hold on;

%color = [0 1 1];
[numRows, numCols] = size(spike_diffs_e_all);

% Initialize a vector to hold the average values
averages_e = zeros(numRows, 1);

% Loop through each row
for row = 1:numRows
    % Extract the values from the current row
    yValues = cell2mat(spike_diffs_e_all(row, :));

    % Calculate the average for the current row
    averages_e(row) = mean(yValues);

    % Initialize colors for the current row
    colors = zeros(length(yValues-1), 3); % m-by-3 matrix for current row

    % Calculate colors for each point in the current row
    for col = 1:length(yValues)
        % Create a gradient from magenta to blue
        colors(col, :) = [(1 - col / (length(yValues) - 1)), 0, (col / (length(yValues) - 1))]; % RGB triplet for magenta to blue
    end

    % Create a scatter plot for the current row
    scatter(conductance_array(row) * ones(size(yValues)), yValues, 100, colors, 'filled');
    hold on; % Keep the current plot
end



    average_array = [];

for i = 1:length(spike_diffs_i_all(:,1))
    average_array = [average_array conductance_array(i)];

end


% Plot the average line
%avgLine_i = plot(average_array, averages_i, 'k-', 'LineWidth', 2, 'DisplayName', 'Mean sequence duration');


% Change font properties
ax = gca; % Get current axes
ax.FontName = 'Helvetica-Narrow'; % Set font name
ax.FontSize = 14; % Set font size

colors = colors(1:end-2,:);

% Create a color bar
clim([1 numCols]); % Set color axis limits
colormap(colors); % Set the colormap to the calculated colors
cb = colorbar; % Create the color bar

% Set the color bar labels to indicate the first and last columns
cb.Ticks = [1, numCols]; % Set ticks at the start and end
cb.TickLabels = {'Earlier sequence', sprintf('Later sequence', numCols)}; % Label the ticks


hold off;
grid on;
xlabel('Conductance [mS/cm^2]');
ylabel('Sequence duration [ms]');
xlim([conductance_array(1)-step_size conductance_array(end)+step_size]); 
title('Sequence durations as function of inhibitory conductance, excitatory population');
xticks(conductance_array(1):step_size:conductance_array(end));

%legend(avgLine_i); % This will display the legend with 'Averages'


%toc