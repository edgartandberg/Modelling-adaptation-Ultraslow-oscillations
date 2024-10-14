function [connectivitymatrix] = plotConnectivityMatrix(spiketrains)

    A = zeros(spiketrains);
    for i = 1:spiketrains
        if i > 1
            A(i, i-1) = 1; % Set the left adjacent element to 1
        end
        if i < spiketrains
            A(i, i+1) = 1; % Set the right adjacent element to 1
        end
    end
    connectivitymatrix = A; 

    heatmap(connectivitymatrix, 'ColorMap', parula, 'GridVisible', 'off');
    title('Connectivity Matrix');
    xlabel('Neuron');
    ylabel('Neuron');
end