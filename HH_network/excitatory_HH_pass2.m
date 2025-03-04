function [v_e, t_e, e_counter, m_e, h_e, n_e, spk_times_e] = excitatory_HH_pass2(v_e, k, t_final,dt,i_ext_e, m_e, h_e, n_e)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


c = 1; % Capacitance
g_k = 10; % Conductance for potassium
g_na = 120; % Conductance for sodium
g_l = 0.1; % Leak conductance
v_k = -70; % Potassium reversal potential
v_na = 55; % Sodium reversal potential
v_l = -65; % Leak reversal potential


spk_times_e = [];

dt05=dt/2;
m_steps=round(t_final/dt);


% potential for excitatory neuron, initial condition
v_e(1)=-70; % potential for excitatory neuron, initial condition



e_counter=0;



v_inc=(g_na*m_e(k)^3*h_e(k)*(v_na-v_e(k))+ g_k*n_e(k)^4*(v_k-v_e(k))+g_l*(v_l-v_e(k))+i_ext_e)/c;
m_inc=alpha_m(v_e(k))*(1-m_e(k))-beta_m(v_e(k))*m_e(k);
h_inc=alpha_h(v_e(k))*(1-h_e(k))-beta_h(v_e(k))*h_e(k);
n_inc=alpha_n(v_e(k))*(1-n_e(k))-beta_n(v_e(k))*n_e(k);
    
v_tmp=v_e(k)+dt05*v_inc;
m_tmp=m_e(k)+dt05*m_inc;
h_tmp=h_e(k)+dt05*h_inc;
n_tmp=n_e(k)+dt05*n_inc;
    
v_inc=(g_na*m_tmp^3*h_tmp*(v_na-v_tmp)+ ...
    g_k*n_tmp^4*(v_k-v_tmp)+g_l*(v_l-v_tmp)+i_ext_e)/c;
m_inc=alpha_m(v_tmp)*(1-m_tmp)-beta_m(v_tmp)*m_tmp;
h_inc=alpha_h(v_tmp)*(1-h_tmp)-beta_h(v_tmp)*h_tmp;
n_inc=alpha_n(v_tmp)*(1-n_tmp)-beta_n(v_tmp)*n_tmp;




if length(v_inc) > 1
    v_e(k+1)=v_e(k)+dt*v_inc(k);
    m_e(k+1)=m_e(k)+dt*m_inc(k);
    h_e(k+1)=h_e(k)+dt*h_inc(k);
    n_e(k+1)=n_e(k)+dt*n_inc(k);

else
    v_e(k+1)=v_e(k)+dt*v_inc;
    m_e(k+1)=m_e(k)+dt*m_inc;
    h_e(k+1)=h_e(k)+dt*h_inc;
    n_e(k+1)=n_e(k)+dt*n_inc;
end


if v_e(k) > 20
    if v_e(k+1) < 20
        e_counter = e_counter+1;
        spk_times_e(e_counter)=k; % make function output this
    end
end




t_e=(0:m_steps)*dt;

end