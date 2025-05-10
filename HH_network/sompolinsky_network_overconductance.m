%% calculate durations as function of conductances

% Load data from Excel file
% filename = 'sequence durations excitatory, simulation over g_ee-kopi.xlsx';
% data = readtable(filename);

%% linear model

% Extract conductance values (Column 1) and averages (Column 5)
conductance = data{:, 1}; % Assuming Column 1 is conductance
conductance = cell2mat(conductance);
averages = data{:, 5};     % Assuming Column 5 is averages

% Fit a linear model
mdl = fitlm(conductance, averages);

% Display the model summary
disp(mdl);

% Interpolate data using the fitted model
% Define new conductance values for interpolation
new_conductance = linspace(min(conductance), max(conductance), 100);
predicted_averages = predict(mdl, new_conductance');

% Plot the results
figure;
scatter(conductance, averages, 'filled'); % Original data points
hold on;
plot(new_conductance, predicted_averages, 'r-', 'LineWidth', 2); % Interpolated line
xlabel('Conductance Values');
ylabel('Averages');
title('Linear Interpolation of Averages');
legend('Data Points', 'Linear Fit');
grid on;

%%

%% plotting excitatory sequence duration as a function of excitatory conductance

% % Load the Excel file
% data = readtable('sequence durations excitatory, simulation over g_ee-kopi.xlsx');

% Extract the first column for x-axis
x = cell2mat(data{:, 1});

% Create a figure
figure;
hold on; % Hold on to plot multiple series

% Get the autumn colormap
colors = autumn(3); % Get 3 colors from the autumn colormap

markerSize = 75; % Increase this value to make points larger


% Plot each column with the autumn colormap
scatter(x, data{:, 2}, markerSize, colors(1, :), 'filled'); % Column 2 in first autumn color
scatter(x, data{:, 3}, markerSize, colors(2, :), 'filled'); % Column 3 in second autumn color
scatter(x, data{:, 4}, markerSize, colors(3, :), 'filled'); % Column 4 in third autumn color

% Add labels and title
xlabel('Excitatory conductance [mS / cm^2]');
ylabel('Sequence Durations [ms]');
title('Excitatory sequence Durations as function of excitatory conductance');
xlim([0.3 1.5]);
ylim([0 300]);


% xlim([0 17]);
% ylim([0 300]);

% Create a color gradient legend
c = colorbar; % Create colorbar
colormap(autumn); % Set the colormap for the colorbar to autumn
% c.Label.String = 'Color Gradient';
c.Ticks = [0, 1]; % Set ticks for the colorbar
c.TickLabels = {'Earlier sequence', 'Later sequence'}; % Label the ticks

hold off; % Release the hold


%% average excitatory sequence duration

% data = readtable('sequence durations excitatory, simulation over g_ee-kopi.xlsx');


% Extract columns
x = cell2mat(data{:, 1}); % Column 1
averages = data{:, 5}; % Column 5 (averages)
% values2 = data{:, 2}; % Column 2
% values3 = data{:, 3}; % Column 3
% values4 = data{:, 4}; % Column 4
% 
% % Calculate variance for each value in column 1
% variance2 = var(values2 - averages); % Variance for column 2
% variance3 = var(values3 - averages); % Variance for column 3
% variance4 = var(values4 - averages); % Variance for column 4

% Create a figure
figure;

% Plot averages
plot(x, averages, 'r-o', 'LineWidth', 2);
xlabel('Excitatory conductance [mS / cm^2]');
ylabel('Average sequence duration [ms]');
title('Average excitatory sequence duration as function of excitatory conductance');
xlim([0.4 1.5]);
ylim([0 300]);
grid on;

%% plotting inhibitory sequence duration as a function of excitatory conductance

% % Load the Excel file
% data = readtable('sequence durations inhibitory, simulation over g_ee-kopi.xlsx');

% Extract the first column for x-axis
x = cell2mat(data{:, 1});

% Create a figure
figure;
hold on; % Hold on to plot multiple series

% Get the autumn colormap
colors = autumn(3); % Get 3 colors from the autumn colormap

markerSize = 75; % Increase this value to make points larger


% Plot each column with the autumn colormap
scatter(x, data{:, 2}, markerSize, colors(1, :), 'filled'); % Column 2 in first autumn color
scatter(x, data{:, 3}, markerSize, colors(2, :), 'filled'); % Column 3 in second autumn color
scatter(x, data{:, 4}, markerSize, colors(3, :), 'filled'); % Column 4 in third autumn color

% Add labels and title
xlabel('Excitatory conductance [mS / cm^2]');
ylabel('Sequence Durations [ms]');
title('Inhibitory sequence Durations as function of excitatory conductance');
xlim([0.3 1.5]);
ylim([0 300]);


% xlim([0 17]);
% ylim([0 300]);

% Create a color gradient legend
c = colorbar; % Create colorbar
colormap(autumn); % Set the colormap for the colorbar to autumn
% c.Label.String = 'Color Gradient';
c.Ticks = [0, 1]; % Set ticks for the colorbar
c.TickLabels = {'Earlier sequence', 'Later sequence'}; % Label the ticks

hold off; % Release the hold

%% average inhibitory sequence duration

% 
% data = readtable('sequence durations inhibitory, simulation over input current.xlsx');


% Extract columns
x = cell2mat(data{:, 1}); % Column 1
averages = data{:, 5}; % Column 5 (averages)
% values2 = data{:, 2}; % Column 2
% values3 = data{:, 3}; % Column 3
% values4 = data{:, 4}; % Column 4
% 
% % Calculate variance for each value in column 1
% variance2 = var(values2 - averages); % Variance for column 2
% variance3 = var(values3 - averages); % Variance for column 3
% variance4 = var(values4 - averages); % Variance for column 4

% Create a figure
figure;

% Plot averages
plot(x, averages, 'b-o', 'LineWidth', 2);
xlabel('Excitatory conductance [mS / cm^2]');
ylabel('Average sequence duration [ms]');
title('Average inhibitory sequence duration as function of excitatory conductance');
xlim([0.4 1.5]);
ylim([0 300]);
grid on;



%%


%  ------   Inhibitory conductance

%% plotting excitatory sequence duration as a function of inhibitory conductance

% % Load the Excel file
% data = readtable('sequence durations excitatory, simulation over g_ii-kopi.xlsx');

% Extract the first column for x-axis
x = cell2mat(data{:, 1});

% Create a figure
figure;
hold on; % Hold on to plot multiple series

% Get the autumn colormap
colors = autumn(3); % Get 3 colors from the autumn colormap

markerSize = 75; % Increase this value to make points larger


% Plot each column with the autumn colormap
scatter(x, data{:, 2}, markerSize, colors(1, :), 'filled'); % Column 2 in first autumn color
scatter(x, data{:, 3}, markerSize, colors(2, :), 'filled'); % Column 3 in second autumn color
scatter(x, data{:, 4}, markerSize, colors(3, :), 'filled'); % Column 4 in third autumn color

% Add labels and title
xlabel('Inhibitory conductance [mS / cm^2]');
ylabel('Sequence Durations [ms]');
title('Excitatory sequence Durations as function of excitatory conductance');
xlim([0 0.12]);
ylim([-50 300]);


% xlim([0 17]);
% ylim([0 300]);

% Create a color gradient legend
c = colorbar; % Create colorbar
colormap(autumn); % Set the colormap for the colorbar to autumn
% c.Label.String = 'Color Gradient';
c.Ticks = [0, 1]; % Set ticks for the colorbar
c.TickLabels = {'Earlier sequence', 'Later sequence'}; % Label the ticks

hold off; % Release the hold


%% average excitatory sequence duration over inhibitory conductance

data = readtable('sequence durations excitatory, simulation over g_ee-kopi.xlsx');


% Extract columns
x = cell2mat(data{:, 1}); % Column 1
averages = data{:, 5}; % Column 5 (averages)
% values2 = data{:, 2}; % Column 2
% values3 = data{:, 3}; % Column 3
% values4 = data{:, 4}; % Column 4
% 
% % Calculate variance for each value in column 1
% variance2 = var(values2 - averages); % Variance for column 2
% variance3 = var(values3 - averages); % Variance for column 3
% variance4 = var(values4 - averages); % Variance for column 4

% Create a figure
figure;

% Plot averages
plot(x, averages, 'r-o', 'LineWidth', 2);
xlabel('Excitatory conductance [mS / cm^2]');
ylabel('Average sequence duration [ms]');
title('Average excitatory sequence duration as function of inhibitory conductance');
xlim([0 0.12]);
ylim([-50 300]);
grid on;

%% plotting inhibitory sequence duration as a function of inhibitory conductance


% % Load the Excel file
% data = readtable('sequence durations inhibitory, simulation over g_ii-kopi.xlsx');

% Extract the first column for x-axis
x = cell2mat(data{:, 1});

% Create a figure
figure;
hold on; % Hold on to plot multiple series

% Get the autumn colormap
colors = autumn(3); % Get 3 colors from the autumn colormap

markerSize = 75; % Increase this value to make points larger


% Plot each column with the autumn colormap
scatter(x, data{:, 2}, markerSize, colors(1, :), 'filled'); % Column 2 in first autumn color
scatter(x, data{:, 3}, markerSize, colors(2, :), 'filled'); % Column 3 in second autumn color
scatter(x, data{:, 4}, markerSize, colors(3, :), 'filled'); % Column 4 in third autumn color

% Add labels and title
xlabel('Inhibitory conductance [mS / cm^2]');
ylabel('Sequence Durations [ms]');
title('Inhibitory sequence Durations as function of inhibitory conductance');
xlim([0 0.12]);
ylim([-50 300]);


% xlim([0 17]);
% ylim([0 300]);

% Create a color gradient legend
c = colorbar; % Create colorbar
colormap(autumn); % Set the colormap for the colorbar to autumn
% c.Label.String = 'Color Gradient';
c.Ticks = [0, 1]; % Set ticks for the colorbar
c.TickLabels = {'Earlier sequence', 'Later sequence'}; % Label the ticks

hold off; % Release the hold

%% average inhibitory sequence duration over inhibitory conductance


% 
% data = readtable('sequence durations inhibitory, simulation over input current.xlsx');
% 

% Extract columns
x = cell2mat(data{:, 1}); % Column 1
averages = data{:, 5}; % Column 5 (averages)
% values2 = data{:, 2}; % Column 2
% values3 = data{:, 3}; % Column 3
% values4 = data{:, 4}; % Column 4
% 
% % Calculate variance for each value in column 1
% variance2 = var(values2 - averages); % Variance for column 2
% variance3 = var(values3 - averages); % Variance for column 3
% variance4 = var(values4 - averages); % Variance for column 4

% Create a figure
figure;

% Plot averages
plot(x, averages, 'b-o', 'LineWidth', 2);
xlabel('Inhibitory conductance [mS / cm^2]');
ylabel('Average sequence duration [ms]');
title('Average inhibitory sequence duration as function of inhibitory conductance');
xlim([0 0.12]);
ylim([-50 300]);
grid on;