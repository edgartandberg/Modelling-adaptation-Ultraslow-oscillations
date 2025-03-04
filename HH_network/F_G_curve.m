% HH network
% conductance-firing rate curve for postsynaptic excitatory neuron

tic


t_final=500; % in ms
dt=0.01;
max_amps = 10; % in mA
e_Hz=zeros(1,max_amps);
i_ext_e = max_amps;
i_ext_i = max_amps;

g_values = (1e-3:1e-3:10e-3); %every step takes approx. 20 s for a simulation over 200ms



% input to first neurons, with external input
[v_e, t_e, e_counter, m, h, n, spk_times_e] = excitatory_HH(t_final,dt,i_ext_e);
[v_i, t_i, i_counter, spk_times_i] = inhibitory_HH(t_final,dt,i_ext_i);



v_e=v_e(1:round(t_final/dt));
v_i=v_i(1:round(t_final/dt));

I_syn_e_all = cell(1,length(g_values));

count = 0;

for g_bar = g_values
   
[I_syn_e, g_e] = synaptic_current_FG(v_e, v_i,t_final,dt, spk_times_e, spk_times_i, g_bar);
count = count+1;

I_syn_e_all(1,count) = {I_syn_e};

end

t_e=t_e(1:round(t_final/dt));
t_i=t_i(1:round(t_final/dt));


%% Connected neurons, taking in post-synaptic current

count=0;
firing_rates = zeros(1, length(g_values));
for k = 1:length(g_values)
    count=count+1;
    disp(count)
    I_syn_e = I_syn_e_all{1, k};
    %disp(I_syn_e)
    [v_e2, t_e, e_counter, m, h, n] = excitatory_HH(t_final,dt,I_syn_e); % excitatory -> excitatory
    Hz_g = e_counter/t_final;
    firing_rates(count)= Hz_g;
end 

[v_i2, t_i, i_counter, m, h, n] = excitatory_HH(t_final,dt,I_syn_i); % inhibitory -> excitatory


%% Plot F_G curve

plot(g_values,firing_rates,'-o','LineWidth',2)
axis([0 1e-2 0 0.11])
title('Firing rate as a function of conductance');
set(gca,'Fontsize',12);
xlabel('g [pA]','Fontsize',20);
ylabel('Firing rate [Hz]','Fontsize',20);




%%
toc
