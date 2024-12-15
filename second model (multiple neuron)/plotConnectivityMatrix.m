function [connectivitymatrix] = plotConnectivityMatrix(spiketrains)
% Takes in spiketrains and spits out connectivity
% ------------------------
% Loops through the spiketrains and makes a n x n matrix, with neurons
% connected to previous one and next on or 'right and left'.
% Also adds noise to weights between neurons, and plots the 
% whole thing as a heatmap.
% ------------------------

% percentage error/noise for weights, change to 0 for weights = 1
error = 0; 

    A = zeros(spiketrains);
    for i = 1:spiketrains
        noise_matrix = 1 - error + error*rand();
        if i > 1
            A(i, i-1) = 1*noise_matrix; 
        end
        if i < spiketrains
            A(i, i+1) = 1*noise_matrix; 
        end
    end
    connectivitymatrix = A; 

    heatmap(connectivitymatrix, 'ColorMap', parula, 'GridVisible', 'off','CellLabelColor','none');
    title('Connectivity Matrix');
    xlabel('Neuron');
    ylabel('Neuron');
end