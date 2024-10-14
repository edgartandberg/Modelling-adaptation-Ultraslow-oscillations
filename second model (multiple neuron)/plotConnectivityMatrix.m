function [connectivitymatrix] = plotConnectivityMatrix(spiketrains)




    A = zeros(spiketrains);
    for i = 1:spiketrains
        noise_matrix = 0.95 + (1 - 0.95) * rand();
        if i > 1
            A(i, i-1) = 1*noise_matrix; % Set to number between 0.5 and 1
        end
        if i < spiketrains
            A(i, i+1) = 1*noise_matrix; % Set to number between 0.5 and 1
        end
    end
    connectivitymatrix = A; 

    heatmap(connectivitymatrix, 'ColorMap', parula, 'GridVisible', 'off');
    title('Connectivity Matrix');
    xlabel('Neuron');
    ylabel('Neuron');
end