clearvars

t_final=200; % in ms
dt=0.01;

i_ext_e = 20; % in micro A
i_ext_i = 20; % in micro A


% [v_e, t_e, e_counter, m, h, n, spiketimes_e] = excitatory_HH(t_final,dt, i_ext_e);
% [v_i, t_i, i_counter, spiketimes_i] = inhibitory_HH(t_final,dt, i_ext_i);

[v_e, t_e, e_counter, m, h, n, spiketimes_e, w_e] = excitatory_HH_adaptation(t_final,dt,i_ext_e);
[v_i, t_i, i_counter, spiketimes_i, w_i] = inhibitory_HH_adaptation(t_final,dt,i_ext_i);



f = figure;

% Set the figure size (width, height)
f.Position = [400, 200, 1200, 600]; % [left, bottom, width, height]

subplot(221)

% First y-axis for v_e
yyaxis left
plot(t_e, v_e, '-k', 'LineWidth', 2);
set(gca, 'FontSize', 12);
title(['Excitatory Neuron - External current = ', num2str(i_ext_e), ' mA']);
xlim([-1 t_final])
ylim([-80 60])
xlabel('t [ms]', 'FontSize', 20);
ylabel('v [mV]', 'FontSize', 20);

ax = gca;
ax.YColor = 'k'; % Change left y-axis label color


% Second y-axis for w_e
yyaxis right
plot(t_e, w_e, 'r'); 
ylabel('w [mA]', 'FontSize', 20); 
ax.YColor = 'red'; % Change right y-axis label color
hold off

subplot(222) % ISI plot for excitatory
plot(diff(spiketimes_e), 'b-*')
xlabel('spike no', 'FontSize',20);
ylabel('time [s]','Fontsize',20);
title('Interspike interval for excitatory neuron')


subplot(223)
plot(t_i,v_i,'-k','Linewidth',2);
title(['Inhibitory Neuron - External current = ', num2str(i_ext_i), ' mA']);
set(gca,'Fontsize',12);
xlim([-1 t_final])
ylim([-80 60])
xlabel('t [ms]','Fontsize',20);
ylabel('v [mV]','Fontsize',20);
ax = gca;


% Second y-axis for w_i
yyaxis right
plot(t_i, w_i, 'r'); 
ylabel('w [mA]', 'FontSize', 20); 
ax.YColor = 'red'; % Change right y-axis label color


subplot(224) % ISI plot for inhibitory
plot(diff(spiketimes_i), 'b-*')
xlabel('spike no', 'FontSize',20);
ylabel('time [s]','Fontsize',20);
title('Interspike interval for inhibitory neuron')
