clc
clearvars

t_final=2000; % in ms
dt=0.01;

i_ext_e = 1.2; % in micro A
i_ext_i = 1.2; % in micro A

Amps = zeros(1,200000);
Amps(20000:180000) = i_ext_e;


ii = 1;



% change to these for model without adaptation
% [v_e, t_e, e_counter, m, h, n, spiketimes_e] = excitatory_HH(t_final,dt,ii, Amps);
% [v_i, t_i, i_counter, spiketimes_i] = inhibitory_HH(t_final,dt, ii, Amps);


% Change to these for model with adaptation current, like in LIF
% [v_e, t_e, e_counter, m_e, h_e, n_e, spiketimes_e, w_e,i_ext_e_all] = excitatory_HH_adaptation(t_final,dt,i_ext_e,t_k,b_k);
% [v_i, t_i, i_counter, spiketimes_i, w_i,i_ext_i_all] = inhibitory_HH_adaptation(t_final,dt,i_ext_i,t_k,b_k);

% Change to these for model with adaptation like in sompolinsky

[v_e, t_e, e_counter, m_e, h_e, n_e, spiketimes_e, i_ext_e_all, I_z, z] = E_Sompolinsky_adaptation(t_final,dt, Amps);
[v_i, t_i, i_counter, m_i, h_i, n_i, spiketimes_i, i_ext_i_all, z_i] = I_Sompolinsky_adaptation(t_final,dt, Amps);



spiketimes_e = spiketimes_e*dt;
spiketimes_i = spiketimes_i*dt;



%% Plots for single neuron spike train


figure(1) 
subplot(2,1,1);
plot(t_e(1:200000),Amps,'k','LineWidth',3)
%title('Input current')
xlabel('Time [ms]', 'FontSize', 14)  
ylabel('Input current [\muA/cm^2]', 'FontSize', 14)  
%xlim([0 t_final])
ylim([-2 i_ext_e*1.5])


subplot(2,1,2);
plot(t_e,v_e, 'k', 'LineWidth',2)
%title('Voltage')
xlabel('Time [ms]', 'FontSize', 14)  
ylabel('V [mV]', 'FontSize', 14)  
%xlim([-0.1*N num_end*1.05])
ylim([-75 70])
%yline(u_th,'--k','Threshold')
hold off



%% Plots with adaptation + input current




f = figure;

% Set the figure size (width, height)
f.Position = [400, 200, 1200, 600]; % [left, bottom, width, height]


subplot(231)
plot(t_e(1:200000),Amps, '-k', 'LineWidth', 2)
ylim([-0.5 3])
xlabel('time [ms]','Fontsize',20);
ylabel('Current [mA]','Fontsize',20);
title('Input Current')

subplot(232)
% First y-axis for v_e
yyaxis left
plot(t_e, v_e, '-r', 'LineWidth', 1);
set(gca, 'FontSize', 12);
%title(['Excitatory Neuron - External current = ', num2str(i_ext_e), ' mA']);
title('Excitatory Neuron');
xlim([-1 t_final])
ylim([-80 60])
xlabel('time [ms]', 'FontSize', 20);
ylabel('Voltage [mV]', 'FontSize', 20);

ax = gca;
ax.YColor = 'k'; % Change left y-axis label color


% Second y-axis for I_z
yyaxis right
plot(t_e, z, '-k', 'LineWidth', 2); 
ylabel('z(t)', 'FontSize', 20); 
ax.YColor = 'black'; % Change right y-axis label color
hold off

subplot(233) % ISI plot for excitatory
plot(diff(spiketimes_e), 'b-*')
xlabel('spike no', 'FontSize',20);
ylabel('time [ms]','Fontsize',20);
title('Interspike interval for excitatory neuron')


subplot(234)
plot(t_i(1:200000),Amps, '-k', 'LineWidth', 2)
ylim([-0.5 3])
xlabel('time [ms]','Fontsize',20);
ylabel('Current [\muA]','Fontsize',20);
title('Input Current')

subplot(235)
plot(t_i,v_i,'-b','Linewidth',1);
%title(['Inhibitory Neuron - External current = ', num2str(i_ext_i), ' mA']);
title('Inhibitory Neuron');
set(gca,'Fontsize',12);
xlim([-1 t_final])
ylim([-80 60])
xlabel('time [ms]','Fontsize',20);
ylabel('Voltage [mV]','Fontsize',20);
ax = gca;



subplot(236) % ISI plot for inhibitory
plot(diff(spiketimes_i), 'b-*')
ylim([60 62])
xlabel('spike no', 'FontSize',20);
ylabel('time [ms]','Fontsize',20);
title('Interspike interval for inhibitory neuron')

%% Plots with adaptation

f = figure;

% Set the figure size (width, height)
f.Position = [400, 200, 1200, 600]; % [left, bottom, width, height]


subplot(221)
% First y-axis for v_e
yyaxis left
plot(t_e, v_e, '-r', 'LineWidth', 2);
set(gca, 'FontSize', 12);
%title(['Excitatory Neuron - External current = ', num2str(i_ext_e), ' mA']);
title('Excitatory Neuron');
xlim([-1 t_final])
ylim([-80 60])
xlabel('t [ms]', 'FontSize', 20);
ylabel('v [mV]', 'FontSize', 20);

ax = gca;
ax.YColor = 'k'; % Change left y-axis label color


% Second y-axis for w_e
yyaxis right
plot(t_e, w_e, '-k', 'LineWidth', 2); 
ylabel('z [mA]', 'FontSize', 20); 
ax.YColor = 'k'; % Change right y-axis label color
hold off

subplot(222) % ISI plot for excitatory
plot(diff(spiketimes_e), 'b-*')
xlabel('spike no', 'FontSize',20);
ylabel('time [ms]','Fontsize',20);
title('Interspike interval for excitatory neuron')




subplot(223)
plot(t_i,v_i,'-b','Linewidth',2);
%title(['Inhibitory Neuron - External current = ', num2str(i_ext_i), ' mA']);
title('Inhibitory Neuron');
set(gca,'Fontsize',12);
xlim([-1 t_final])
ylim([-80 60])
xlabel('t [ms]','Fontsize',20);
ylabel('v [mV]','Fontsize',20);
ax = gca;


% Second y-axis for w_i
yyaxis right
plot(t_i, w_i, '-k', 'LineWidth', 2); 
ylabel('w [mA]', 'FontSize', 20); 
ax.YColor = 'k'; % Change right y-axis label color


subplot(224) % ISI plot for inhibitory
plot(diff(spiketimes_i), 'b-*')
xlabel('spike no', 'FontSize',20);
ylabel('time [ms]','Fontsize',20);
title('Interspike interval for inhibitory neuron')