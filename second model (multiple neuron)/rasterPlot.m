function rasterPlot(spikeCountsCell, timeStepS)


% color = lines(length(spikeCountsCell)); % Generates distinct colors. 
color = 'k' ; % switch to this line for single color


% Loop through each cell in the cell array
for cellIndex = 1:length(spikeCountsCell)
    spikes = spikeCountsCell{cellIndex}; % Get the spike train matrix from the current cell
    trains = size(spikes, 1); 
    ticMargin = 0.01;                                      % gap between spike trains (full scale is 1)
    ticHeight = (1.0 - (trains + 1) * ticMargin) / trains;

    for train = 1:trains
        spikeTimes = find(spikes(train, :) == 1);
        yOffset = ticMargin + (train - 1) * (ticMargin + ticHeight) + (cellIndex - 1) * (1 + ticMargin);
        for i = 1:length(spikeTimes)
            line([spikeTimes(i) * timeStepS, spikeTimes(i) * timeStepS], ...
                 [yOffset, yOffset + ticHeight], 'Color', color); % remove '(cellIndex, :)' for single color
        end                                                       % add '(cellIndex, :)' for varying colors    
    end
end

% Set axis labels and title
xlabel('Time (s)');
ylabel('Neuron');
title('Raster plot');
hold off;

end