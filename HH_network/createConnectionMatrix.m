


function connectivity = createConnectionMatrix(n)
    % Validate input
    if mod(n, 3) ~= 0
        error('n must be a multiple of 3.');
    end

    % Initialize the connection matrix
    connectionMatrix = zeros(n, n/3);

    % Calculate the number of connections
    if mod(n, 2) == 0
        leftConnections = n / 6;  % Half of n/3 for left
        rightConnections = n / 6; % Half of n/3 for right
    else
        leftConnections = floor(n / 6);  % Floor for left
        rightConnections = ceil(n / 6);  % Ceiling for right
    end

    for i = 1:n
        % Connect to left neurons
        for j = 1:leftConnections
            connectionIndex = mod(i - j - 1, n) + 1 ; % Left connection
            connectionMatrix(i, j) = connectionIndex;
        end

        % % Connect to right neurons
        % for j = 1:rightConnections
        %     connectionIndex = mod(i + j - 1, n) + 1 ; % Right connection
        %     connectionMatrix(i, leftConnections + j) = connectionIndex;
        % end
    end



    for i = 1:n
        for j = 1:size(connectionMatrix, 2) % Loop over all columns
           connIndex = connectionMatrix(i, j); % Get the connection index
            if connIndex > 0 && connIndex <= n % Check for valid index
                connectivity(i, connIndex) = 1; % Set the connectivity for the ith row
            end
        end
    end

end



