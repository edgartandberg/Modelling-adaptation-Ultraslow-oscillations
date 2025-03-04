function [spike_times_e, spike_times_i] = rasterPlot_HH(spiketrains_e_all,spiketrains_i_all)

timeStepS = 1;

% color = lines(length(spikeCountsCell)); % switxh to this line for distinct colors. 
color = 'k' ; % switch to this line for single color


% Loop through each cell in the cell array
% from the 'Poisson tutorial' :
% https://www.hms.harvard.edu/bss/neuro/bornlab/nb204/statistics/poissonTutorial.txt

% Function to detect spikes

% Initialize new cell arrays to hold spike times
spike_times_e = cell(size(spiketrains_e_all));
spike_times_i = cell(size(spiketrains_i_all));

% Function to detect spikes
function spike_times = detect_spikes(spiketrain)
    spike_times = []; % Initialize spike times array
    for t = 1:length(spiketrain)-1
        if spiketrain(t) < 20 && spiketrain(t+1) > 20
            spike_times = [spike_times, t]; % Record the time of the spike
        end
    end
end

% Process excitatory spiketrains
for neuron = 1:numel(spiketrains_e_all)
    spike_times_e{neuron} = detect_spikes(spiketrains_e_all{neuron});
end

% Process inhibitory spiketrains
for neuron = 1:numel(spiketrains_i_all)
    spike_times_i{neuron} = detect_spikes(spiketrains_i_all{neuron});
end


% for excitatory raster plot
pre_array = zeros(1,length(spiketrains_e_all{1}));


for cellIndex = 1:length(spike_times_e)
    
    
    spikes = spike_times_e{cellIndex}; % Get the spike train matrix from the current cell
    trains = size(spikes, 1); 

    % Convert spike times to indices for pre_array
    for i = 1:length(spikes)
        index = round(spikes(i) * 0.01); % Convert ms to index (assuming 1 ms corresponds to index 100)
        pre_array(index) = 1; % Set the corresponding index in pre_array to 1
    end
    
    spikes = pre_array; % pre-assigning empty array
    ticMargin = 0.01;                                      % gap between spike trains (full scale is 1)
    ticHeight = (0.7 - (trains + 1) * ticMargin) / trains;

    for train = 1:trains

        spikeTimes = find(spikes(train, :) == 1);
        yOffset = ticMargin + (train - 1) * (ticMargin + ticHeight) + (cellIndex - 1) * (1 + ticMargin);
        for i = 1:length(spikeTimes)
            line([spikeTimes(i) * timeStepS, spikeTimes(i) * timeStepS], ...
                 [yOffset, yOffset + ticHeight], 'Color', color, 'linewidth', 2); % remove '(cellIndex, :)' for single color
        end                                                       % add '(cellIndex, :)' for varying colors    
    end
end

% Set axis labels and title
%figure()
xlabel('Time (ms)');
ylabel('Neuron');
title('Raster plot');
%axis([0, length(spikes)/100, 0, numel(spiketrains_e_all) ]);
%hold on;

% % for inhibitory raster plot
% pre_array = zeros(1,length(spiketrains_i_all{1}));
% 
% 
% for cellIndex = 1:length(spike_times_i)
% 
% 
%     spikes = spike_times_i{cellIndex}; % Get the spike train matrix from the current cell
%     trains = size(spikes, 1); 
% 
%     % Convert spike times to indices for pre_array
%     for i = 1:length(spikes)
%         index = round(spikes(i) * 0.01); % Convert ms to index (assuming 1 ms corresponds to index 100)
%         pre_array(index) = 1; % Set the corresponding index in pre_array to 1
%     end
% 
%     spikes = pre_array; % pre-assigning empty array
%     ticMargin = 0.01;                                      % gap between spike trains (full scale is 1)
%     ticHeight = (0.7 - (trains + 1) * ticMargin) / trains;
% 
%     for train = 1:trains
% 
%         spikeTimes = find(spikes(train, :) == 1);
%         yOffset = ticMargin + (train - 1) * (ticMargin + ticHeight) + (cellIndex - 1) * (1 + ticMargin);
%         for i = 1:length(spikeTimes)
%             line([spikeTimes(i) * timeStepS, spikeTimes(i) * timeStepS], ...
%                  [yOffset, yOffset + ticHeight], 'Color', 'b', 'linewidth', 2); % remove '(cellIndex, :)' for single color
%         end                                                       % add '(cellIndex, :)' for varying colors    
%     end
% end




% Set axis labels and title

xlabel('Time (ms)');
ylabel('Neuron');
title('Raster plot');
%axis([0, length(spikes)/100, 0, numel(spiketrains_e_all) ]);
hold on;


end