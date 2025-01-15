% single hodgkin-huxley model

t_final=2000; % in ms
dt=0.01;

max_amps = 20; % in mA

e_Hz=zeros(1,max_amps);
i_Hz=zeros(1,max_amps);

for ii= 1:max_amps

    i_ext_e = ii; % injected current to excitatory neurons
    i_ext_i = ii; % injected current to inhibitory neurons

    [v_e, t_e, e_counter] = excitatory_HH(t_final,dt, i_ext_e);
    [v_i, t_i, i_counter] = inhibitory_HH(t_final,dt, i_ext_i);

    e_Hz(ii) = e_counter/(t_final*0.001);
    i_Hz(ii) = i_counter/(t_final*0.001);

end


figure(1)
plot(e_Hz)


% figure(2)
% subplot(211)
% plot(t_e,v_e,'-k','Linewidth',2);
% title(['excitatory, i = ', num2str(i_ext_e), '    firing rate = ', num2str(e_spikerate)]);
% set(gca,'Fontsize',12);
% xlabel('t [ms]','Fontsize',20);
% ylabel('v [mV]','Fontsize',20);
% 
% 
% subplot(212)
% plot(t_i,v_i,'-k','Linewidth',2);
% title(['inhibitory, i = ', num2str(i_ext_i), '    firing rate = ', num2str(i_spikerate)]);
% set(gca,'Fontsize',12);
% xlabel('t [ms]','Fontsize',20);
% ylabel('v [mV]','Fontsize',20);