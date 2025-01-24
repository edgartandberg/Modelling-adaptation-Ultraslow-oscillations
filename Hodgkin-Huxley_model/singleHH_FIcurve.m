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
clearvars

t_final=2000; % in ms
dt=0.01;

max_amps = 10; % in micro A
step_size_amps = 0.1; % step size for I, in micro A

Amps = 0:step_size_amps:max_amps;

% e_Hz=zeros(1,length(Amps));
% i_Hz=zeros(1,length(Amps));

for ii = 1:length(Amps)

    i_ext_e = Amps(ii); % injected current to excitatory neurons
    i_ext_i = Amps(ii); % injected current to inhibitory neurons
   

    [v_e, t_e, e_counter, m, h, n] = excitatory_HH(t_final,dt, i_ext_e);
    [v_i, t_i, i_counter] = inhibitory_HH(t_final,dt, i_ext_i);

    e_Hz(ii) = e_counter/(t_final*0.001);
    i_Hz(ii) = i_counter/(t_final*0.001);

end

lineColor = [0, 0.4470, 0.3410]; % Emerald green

% Create a new figure
figure;

% First subplot for excitatory neuron
subplot(1, 2, 1); % 1 row, 2 columns, first subplot
plot(Amps, e_Hz, '-*', 'LineWidth', 4, 'Color', lineColor);
title('Excitatory neuron');
xlim([0 max_amps])
xlabel('I [mA / cm^2]', 'FontSize', 15);
ylabel('Firing rate [Hz]', 'FontSize', 15);

% Second subplot for inhibitory neuron
subplot(1, 2, 2); % 1 row, 2 columns, second subplot
plot(Amps, i_Hz, '-*', 'LineWidth', 4, 'Color', lineColor);
title('Inhibitory neuron');
xlim([0 max_amps])
xlabel('I [mA / cm^2]', 'FontSize', 15);
ylabel('Firing rate [Hz]', 'FontSize', 15);

% Adjust layout
sgtitle('Gain Functions of Hodgkin-Huxley Model'); % Overall title for the figure



figure(2)
subplot(211)
plot(t_e,v_e,'-k','Linewidth',2);
%title(['excitatory, i = ', num2str(i_ext_e), '    firing rate = ', num2str(e_Hz)]);
set(gca,'Fontsize',12);
xlim([-1 100])
ylim([-80 50])
xlabel('t [ms]','Fontsize',20);
ylabel('v [mV]','Fontsize',20);


subplot(212)
plot(t_i,v_i,'-k','Linewidth',2);
%title(['inhibitory, i = ', num2str(i_ext_i), '    firing rate = ', num2str(i_Hz)]);
set(gca,'Fontsize',12);
xlim([-1 100])
ylim([-70 50])
xlabel('t [ms]','Fontsize',20);
ylabel('v [mV]','Fontsize',20);

%toc

%% Plot for gating variables for excitatory neuron


figure();

% Plot each variable with different colors
yyaxis right
plot(t_e,v_e,'-k','Linewidth',4,'DisplayName', 'v');
%title(['excitatory, i = ', num2str(i_ext_e), '    firing rate = ', num2str(e_Hz)]);
set(gca,'Fontsize',12);
hold on; % Hold on to add more plots
yyaxis left
plot(t_e, m,  'r--', 'DisplayName', 'm'); % Red color for m
plot(t_e, n, 'g--', 'DisplayName', 'n'); % Green color for n
plot(t_e, h,  'b--', 'DisplayName', 'h'); % Blue color for h

xlim([0 30]);

title('Plot of gating variables with membrane potential');
xlabel('Time');
yyaxis left
ylabel('x_0(t)');
yyaxis right
ylabel('v(t)');

legend('show'); % Display the legend
hold off;

%% Plot for sodium and potassium current

figure();

% Plot each variable with different colors
yyaxis right
plot(t_e,v_e,'-k','Linewidth',4,'DisplayName', 'v');
%title(['excitatory, i = ', num2str(i_ext_e), '    firing rate = ', num2str(e_Hz)]);
set(gca,'Fontsize',12);
hold on; % Hold on to add more plots

% sodium current = (g_na*m(k)^3*h(k)*(v_na-v_e(k))
% potassium current = g_k*n(k)^4*(v_k-v_e(k))

yyaxis left
plot(t_e, m,  'r--', 'DisplayName', 'm'); % Red color for m
plot(t_e, n, 'g--', 'DisplayName', 'n'); % Green color for n
plot(t_e, h,  'b--', 'DisplayName', 'h'); % Blue color for h

xlim([0 30]);

title('Plot of potassium and sodium current');
xlabel('Time');
yyaxis left
ylabel('x_0(t)');
yyaxis right
ylabel('v(t)');

legend('show'); % Display the legend
hold off;