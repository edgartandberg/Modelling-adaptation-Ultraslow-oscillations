% Function to detect spikes
function spike_times = detect_spikes(spiketrain)
    spike_times = []; % Initialize spike times array
    for t = 1:length(spiketrain)-1
        if spiketrain(t) < 20 && spiketrain(t+1) > 20
            spike_times = [spike_times, t]; % Record the time of the spike
        end
    end
end