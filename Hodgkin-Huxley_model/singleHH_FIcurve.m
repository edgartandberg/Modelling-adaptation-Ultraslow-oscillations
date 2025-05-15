% single hodgkin-huxley model

% for 
% t_final=20000; % in ms
% dt=0.001;
% simulation time is about 910 seconds

% for
% t_final=2000; % in ms
% dt=0.01;
% simulation takes around 9 seconds

%tic
clc
clearvars

t_final=2000; % in ms
dt=0.01;

m_steps = round(t_final/dt);

max_amps = 8; % in micro A
step_size_amps = 0.01; % step size for I, in micro A

Amps = 0:step_size_amps:max_amps;

i_ext_e_all = zeros(length(Amps), m_steps);
i_ext_i_all = zeros(length(Amps), m_steps);


% e_Hz=zeros(1,length(Amps));
% i_Hz=zeros(1,length(Amps));

for ii = 1:length(Amps)


    [v_e, t_e, e_counter, m_e, h_e, n_e, k, i_ext_e_all] = excitatory_HH(t_final,dt, ii, Amps);
    [v_i, t_i, i_counter, spiketimes_i, w_i,i_ext_i_all] = inhibitory_HH(t_final,dt, ii, Amps);


    e_Hz(ii) = e_counter/(t_final*0.001);
    i_Hz(ii) = i_counter/(t_final*0.001);

end

lineColor = [0, 0.4470, 0.3410]; % Emerald green

%% Plot for single gain function

plot(Amps(1:748), e_Hz, '-*', 'LineWidth', 2, 'Color','k');
title('Excitatory neuron');
xlim([0 max_amps])
xlabel('I [mA / cm^2]', 'FontSize', 15);
ylabel('Firing rate [Hz]', 'FontSize', 15);


%% plot for gain functions E and I

% Create a new figure
figure;

% First subplot for excitatory neuron
subplot(1, 2, 1); % 1 row, 2 columns, first subplot
plot(Amps, e_Hz, '-*', 'LineWidth', 4, 'Color','r');
title('Excitatory neuron');
xlim([0 max_amps])
xlabel('I [mA / cm^2]', 'FontSize', 15);
ylabel('Firing rate [Hz]', 'FontSize', 15);

% Second subplot for inhibitory neuron
subplot(1, 2, 2); % 1 row, 2 columns, second subplot
plot(Amps, i_Hz, '-*', 'LineWidth', 4, 'Color', 'b');
title('Inhibitory neuron');
xlim([0 max_amps])
xlabel('I [mA / cm^2]', 'FontSize', 15);
ylabel('Firing rate [Hz]', 'FontSize', 15);

% Adjust layout
sgtitle('Gain Functions for Hodgkin-Huxley Model'); % Overall title for the figure

%% plot of spike trains

figure(2)
subplot(211)
yyaxis left
plot(t_e,v_e,'-r','Linewidth',2);
set(gca, 'YColor', 'k'); % Set left y-axis color to black
%title(['excitatory, i = ', num2str(i_ext_e), '    firing rate = ', num2str(e_Hz)]);
set(gca,'Fontsize',12);
xlim([0 500])
ylim([-70 70])
xlabel('t [ms]','Fontsize',20);
ylabel('v [mV]','Fontsize',20);
%yline(i_ext_e,'r','LineWidth',1)
yyaxis right
plot(i_ext_e_all(end,:),'k','Linewidth',2)
ylim([0 10])
set(gca, 'YColor', 'k'); % Set right y-axis color to red
ylabel('Current [mA]','Fontsize',20);



subplot(212)
yyaxis left
plot(t_i,v_i,'-b','Linewidth',2);
set(gca, 'YColor', 'k'); % Set left y-axis color to black
%title(['inhibitory, i = ', num2str(i_ext_i), '    firing rate = ', num2str(i_Hz)]);
set(gca,'Fontsize',12);
xlim([0 500])
ylim([-70 70])
xlabel('t [ms]','Fontsize',20);
ylabel('v [mV]','Fontsize',20);
%yline(i_ext_i,'r','LineWidth',1)
yyaxis right
plot(i_ext_i_all(end,:),'k','Linewidth',2)
ylim([0 10])
set(gca, 'YColor', 'k'); % Set right y-axis color to red
ylabel('Current [mA]','Fontsize',20);


%toc

%% Plot for gating variables for excitatory neuron


figure();

% Plot each variable with different colors
yyaxis right
plot(t_e,v_e,'-r','Linewidth',2,'DisplayName', 'v');
%title(['excitatory, i = ', num2str(i_ext_e), '    firing rate = ', num2str(e_Hz)]);
set(gca,'Fontsize',12);
hold on; % Hold on to add more plots
yyaxis left
plot(t_e, m_e, 'Color', "#0072BD", 'LineStyle', '--', 'LineWidth', 1, 'DisplayName', 'm'); % Magenta
plot(t_e, n_e, 'Color', "g", 'LineStyle', '--', 'LineWidth', 1, 'DisplayName', 'n'); % Green
plot(t_e, h_e, 'Color', "#EDB120", 'LineStyle', '--', 'LineWidth', 1,'DisplayName', 'h'); % Blue


xlim([0 100]);

title('Gating variables for excitatory neuron');
xlabel('Time');
yyaxis left
ylabel('x(t)');
set(gca, 'YColor', 'k'); % Set right y-axis color to red

yyaxis right
ylabel('Voltage [mV]');
set(gca, 'YColor', 'r'); % Set right y-axis color to red


legend('show'); % Display the legend
hold off;



%% Plot for gating variables for inhibitory neuron
figure();

% Plot each variable with different colors
yyaxis right
plot(t_i,v_i,'-b','Linewidth',1,'DisplayName', 'v');
%title(['excitatory, i = ', num2str(i_ext_e), '    firing rate = ', num2str(e_Hz)]);
set(gca,'Fontsize',12);
hold on; % Hold on to add more plots

% sodium current = (g_na*m(k)^3*h(k)*(v_na-v_e(k))
% potassium current = g_k*n(k)^4*(v_k-v_e(k))

yyaxis left
plot(t_e, m_e, 'Color', "y", 'LineStyle', '--', 'LineWidth', 1, 'DisplayName', 'm'); % Magenta
plot(t_e, n_e, 'Color', "g", 'LineStyle', '--', 'LineWidth', 1, 'DisplayName', 'n'); % Green
plot(t_e, h_e, 'Color', "r", 'LineStyle', '--', 'LineWidth', 1,'DisplayName', 'h'); % Blue

xlim([0 100]);

title('Gating variables for inhibitory neuron');
xlabel('Time');
yyaxis left
ylabel('x(t)');
set(gca, 'YColor', 'k'); % Set right y-axis color to red

yyaxis right
ylabel('Voltage [mV]');
set(gca, 'YColor', 'b'); % Set right y-axis color to black


legend('show'); % Display the legend
hold off;