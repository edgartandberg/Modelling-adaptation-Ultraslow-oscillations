clc
clearvars

t_final=500; % in ms
dt=0.01;

i_ext_e = 10.0; % in micro A
i_ext_i = 10.0; % in micro A

max_tau = 150;
max_bk = 0.015;

t_k_values = 10:10:max_tau; 
b_k_values = 0.001:0.001:max_bk; % b_k from 0 to 60 in steps of 5

e_hz = zeros(length(b_k_values),length(t_k_values));
i_hz = zeros(length(b_k_values),length(t_k_values));



for i = 1:length(t_k_values)
    t_k = t_k_values(i);
    for j = 1:length(b_k_values)
        b_k = b_k_values(j);
[v_e, t_e, e_counter, m, h, n, spiketimes_e, w_e,i_ext_e_all] = excitatory_HH_adaptation(t_final,dt,i_ext_e,t_k,b_k);
[v_i, t_i, i_counter, spiketimes_i, w_i, i_ext_i_all] = inhibitory_HH_adaptation(t_final,dt,i_ext_i,t_k,b_k);

e_hz(j,i) = length(spiketimes_e) / t_final;
i_hz(j,i) = length(spiketimes_i) / t_final;
    end
end



spiketimes_e = spiketimes_e*dt;
spiketimes_i = spiketimes_i*dt;


%% Plot heatmap for E neurons

Xlabels = t_k_values;
Xlabels(mod(Xlabels,50) ~= 0) = " ";

Ylabels = b_k_values;
Ylabels(mod(Ylabels,0.004) ~= 0) = " ";

figure()
h = heatmap(t_k_values,b_k_values,e_hz,'CellLabelColor','none');
%caxis([-0.0 2.8]);
h.XDisplayLabels = Xlabels;
h.YDisplayLabels = Ylabels;

h.Title = 'Firing rate for varying time constant and adaptation gain, excitatory neurons';
h.XLabel = 'time constant';
h.YLabel = 'adaptation gain';
h.Colormap = autumn;

set(gca, 'YDir','reverse')


%clim([0 2.5])


%% Plot heatmap for I neurons

Xlabels = t_k_values;
Xlabels(mod(Xlabels,50) ~= 0) = " ";

Ylabels = b_k_values;
Ylabels(mod(Ylabels,0.004) ~= 0) = " ";

figure()
h = heatmap(t_k_values,b_k_values,i_hz,'CellLabelColor','none');
%caxis([-0.0 2.8]);
h.XDisplayLabels = Xlabels;
h.YDisplayLabels = Ylabels;

h.Title = 'Firing rate for varying time constant and adaptation gain, inhibitory neurons';
h.XLabel = 'time constant';
h.YLabel = 'adaptation gain';
h.Colormap = autumn;

h.YDir = 'reverse'; % This line flips the y-axis

%clim([0 2.5])