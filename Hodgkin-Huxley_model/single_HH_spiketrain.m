clearvars

t_final=2000; % in ms
dt=0.01;

i_ext_e = 3; % in micro A
i_ext_i = 3; % in micro A


[v_e, t_e, e_counter, m, h, n] = excitatory_HH(t_final,dt, i_ext_e);
[v_i, t_i, i_counter] = inhibitory_HH(t_final,dt, i_ext_i);


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