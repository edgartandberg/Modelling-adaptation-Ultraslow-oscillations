% HH network
% Network model with basic HH single-neuron model

tic

t_final=200; % in ms
dt=0.01;
max_amps = 10; % in mA
e_Hz=zeros(1,max_amps);
i_ext_e = max_amps;
i_ext_i = max_amps;

size_network = 6; % amount of both excitatory and inhibitory neurons, size of network

%% input to first neurons, with external input


spiketrains_e_all = cell(1, size_network);
spiketrains_i_all = cell(1, size_network);

for j = 1:size_network;

    [v_e, t_e, e_counter, m, h, n, spk_times_e] = excitatory_HH(t_final,dt,i_ext_e);
    [v_i, t_i, i_counter, spk_times_i] = inhibitory_HH(t_final,dt,i_ext_i);

    spiketrains_e_all{j} = v_e;
    spiketrains_i_all{j} = v_i;

end

excitatory_neurons = (1:1:size_network);
inhibitory_neurons = (1:1:size_network);

r = randi([1 size_network],1,round(size_network/3));


%spiketrains_e_all{r} 

%%


e_e_synapses = cell(1,length(r));

for jj = r
    spiketrains_e_all{jj} = v_e;
    [I_syn_e, I_syn_i, g_e] = synaptic_current(v_e, v_i,t_final,dt, spk_times_e, spk_times_i);
    
    e_e_synapses{jj} = I_syn_e;
end

nonEmptyIndices = find(~cellfun(@isempty, e_e_synapses)); %find indices of excitatory that receive e synaptic input


%%


% v_e=v_e(1:round(t_final/dt));
% v_i=v_i(1:round(t_final/dt));
% 
% [I_syn_e, I_syn_i, g_e] = synaptic_current(v_e, v_i,t_final,dt, spk_times_e, spk_times_i);
% 
% t_e=t_e(1:round(t_final/dt));
% t_i=t_i(1:round(t_final/dt));

% figure()
% plot(t_i, I_syn_e, 'DisplayName', 'excitatory', 'LineWidth', 2); % Thicker line for excitatory current
% hold on
% plot(t_i, I_syn_i, 'DisplayName', 'inhibitory', 'LineWidth', 2); % Thicker line for inhibitory current
% axis([0 100 -15 10]); % Axis limitset var 
% legend show; % Display the legend
% title('Plot of Excitatory and Inhibitory Post-Synaptic Current'); % Add title

%% Connected neurons, taking in post-synaptic current

connected_e = spiketrains_e_all(nonEmptyIndices); % voltage of e neurons, that are to receive excitatory input


for ii = 1:connected_e
    connected_e{ii} = v_e;
    [v_e_after_syn, t_e, e_counter, m, h, n, spk_times_e] = excitatory_HH_pass2(v_e, t_final,dt,i_ext_e);

end



[v_i2, t_i, i_counter, m, h, n] = excitatory_HH(t_final,dt,I_syn_i); % inhibitory -> excitatory

%[v_i2, t_i, i_counter] = inhibitory_HH(t_final,dt,I_syn_i); % inhibitory
% -> inhibitory

%% Plot postsynaptic membrane potentials

t_e=t_e(1:round(t_final/dt));
t_i=t_i(1:round(t_final/dt));


figure(2)
subplot(321)
plot(t_e,v_e,'-k','LineWidth',2)
axis([0 200 -75 60])
title('Synapse 1: presynaptic excitatory neuron');
set(gca,'Fontsize',12);
xlabel('t [ms]','Fontsize',20);
ylabel('v [mV]','Fontsize',20);

subplot(322)
plot(t_i,v_i,'-k','LineWidth',2)
axis([0 200 -75 60])
title('Synapse 2: presynaptic inhibitory neuron');
set(gca,'Fontsize',12);
xlabel('t [ms]','Fontsize',20);
ylabel('v [mV]','Fontsize',20);

subplot(323)
plot(t_e, I_syn_e,'-r','LineWidth',2)
axis([0 200 -5 5])
title('input current from presynaptic E to postsynaptic E neuron')
set(gca,'Fontsize',12);
xlabel('t [ms]','Fontsize',20);
ylabel('Current [mA]','Fontsize',20);

subplot(324)
plot(t_i, I_syn_i,'-r','LineWidth',2)
axis([0 200 -5 5])
title('input current from presynaptic I to postsynaptic E neuron')
set(gca,'Fontsize',12);
xlabel('t [ms]','Fontsize',20);
ylabel('Current [mA]','Fontsize',20);


subplot(325)
plot(t_e,v_e2,'-k','Linewidth',2);
axis([0 200 -75 60])
title('Synapse 1: postsynaptic excitatory neuron');
set(gca,'Fontsize',12);
xlabel('t [ms]','Fontsize',20);
ylabel('v [mV]','Fontsize',20);


subplot(326)
plot(t_i,v_i2,'-k','Linewidth',2);
axis([0 200 -75 60])
title('Synapse 2: postsynaptic excitatory neuron');
set(gca,'Fontsize',12);
xlabel('t [ms]','Fontsize',20);
ylabel('v [mV]','Fontsize',20);


%%
toc
