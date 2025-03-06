function [v_i, t_i, i_counter, m_i, h_i, n_i, spk_times_i] = inhibitory_HH_pass2(v_i, k, t_final,dt,i_ext_i, m_i, h_i, n_i)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here


c = 1; % Capacitance
g_k = 20; % Conductance for potassium
g_na = 120; % Conductance for sodium
g_l = 0.1; % Leak conductance
v_k = -70; % Potassium reversal potential
v_na = 55; % Sodium reversal potential
v_l = -65; % Leak reversal potential


spk_times_i = [];

dt05=dt/2;
m_steps=round(t_final/dt);



v_i(1)=-70; % potential for inhibitory neuron, initial condition




i_counter=0;


    
v_inc=(g_na*m_i(k-1)^3*h_i(k-1)*(v_na-v_i(k-1))+ ...
    g_k*n_i(k-1)^4*(v_k-v_i(k-1))+g_l*(v_l-v_i(k-1))+i_ext_i)/c;
m_inc=alpha_m(v_i(k-1))*(1-m_i(k-1))-beta_m(v_i(k-1))*m_i(k-1);
h_inc=alpha_h(v_i(k-1))*(1-h_i(k-1))-beta_h(v_i(k-1))*h_i(k-1);
n_inc=alpha_n(v_i(k-1))*(1-n_i(k-1))-beta_n(v_i(k-1))*n_i(k-1);

v_tmp=v_i(k-1)+dt05*v_inc;
m_tmp=m_i(k-1)+dt05*m_inc;
h_tmp=h_i(k-1)+dt05*h_inc;
n_tmp=n_i(k-1)+dt05*n_inc;

v_inc=(g_na*m_tmp^3*h_tmp*(v_na-v_tmp)+ ...
    g_k*n_tmp^4*(v_k-v_tmp)+g_l*(v_l-v_tmp)+i_ext_i)/c;
m_inc=alpha_m(v_tmp)*(1-m_tmp)-beta_m(v_tmp)*m_tmp;
h_inc=alpha_h(v_tmp)*(1-h_tmp)-beta_h(v_tmp)*h_tmp;
n_inc=alpha_n(v_tmp)*(1-n_tmp)-beta_n(v_tmp)*n_tmp;



if length(v_inc) > 1
    v_i(k)=v_i(k-1)+dt*v_inc(k-1);
    m_i(k)=m_i(k-1)+dt*m_inc(k-1);
    h_i(k)=h_i(k-1)+dt*h_inc(k-1);
    n_i(k)=n_i(k-1)+dt*n_inc(k-1);

else
    v_i(k)=v_i(k-1)+dt*v_inc;
    m_i(k)=m_i(k-1)+dt*m_inc;
    h_i(k)=h_i(k-1)+dt*h_inc;
    n_i(k)=n_i(k-1)+dt*n_inc;
end

if v_i(k-1) > 20
    if v_i(k) < 20
        i_counter = i_counter+1;
        spk_times_i(i_counter)=k;
    end
end




t_i=(0:m_steps)*dt;


end