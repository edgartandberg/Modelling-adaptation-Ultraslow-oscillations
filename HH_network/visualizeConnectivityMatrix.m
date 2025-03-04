function visualizeConnectivityMatrix(connectivity)


    % Visualize  connectivity matrix
    %figure();
    %imagesc(connectivity);

    % Define a custom colormap: black for 0 (not connected), dark blue for 1 (connected)
    colormap([0 0 0; 0 0 0.5]); % White for 0, dark blue for 1
    colorbar('off'); % Turn off the color scale
    title('Connectivity Matrix');
    xlabel('Neuron no');
    ylabel('Neuron no');
    axis square; % Make axes equal

    % Set ticks and labels for both axes
    n = size(connectivity, 1); % Get the number of neurons

    set(gca, 'XTick', [], 'YTick', []); % Remove tick marks

    % set(gca, 'XTick', 1:n, 'YTick', 1:n); % Set ticks for all neurons
    % set(gca, 'XTickLabel', 1:n, 'YTickLabel', 1:n); % Label ticks with neuron indices
end