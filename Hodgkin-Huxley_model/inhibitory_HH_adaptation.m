function [v_i, t_i, i_counter, spiketimes_i, w_i] = inhibitory_HH_adaptation(t_final,dt,i_ext_i)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

dt05=dt/2;


c = 1; % Capacitance
g_k = 20; % Conductance for potassium
g_na = 120; % Conductance for sodium
g_l = 0.1; % Leak conductance
v_k = -70; % Potassium reversal potential
v_na = 55; % Sodium reversal potential
v_l = -65; % Leak reversal potential


t_k = 40; % time constant for adaptation current, indexing over k neurons
a_k = 0.001; % adaptation coupling
b_k = 0.01; % adaptation gain




m_steps=round(t_final/dt);
w = zeros(1,m_steps);
w(1)=0;  % initial conditions adaptation current


v_i(1)=-70; % potential for inhibitory neuron, initial condition

% gating variables, dimensionless
m(1)=alpha_m(v_i(1))/(alpha_m(v_i(1))+beta_m(v_i(1)));
h(1)=0.7; 
n(1)=0.6; 

i_counter=0;

for k=1:m_steps
    
    v_inc=(g_na*m(k)^3*h(k)*(v_na-v_i(k))+ ...
        g_k*n(k)^4*(v_k-v_i(k))+g_l*(v_l-v_i(k))+i_ext_i)/c;
    m_inc=alpha_m(v_i(k))*(1-m(k))-beta_m(v_i(k))*m(k);
    h_inc=alpha_h(v_i(k))*(1-h(k))-beta_h(v_i(k))*h(k);
    n_inc=alpha_n(v_i(k))*(1-n(k))-beta_n(v_i(k))*n(k);
    
    v_tmp=v_i(k)+dt05*v_inc;
    m_tmp=m(k)+dt05*m_inc;
    h_tmp=h(k)+dt05*h_inc;
    n_tmp=n(k)+dt05*n_inc;
    
    v_inc=(g_na*m_tmp^3*h_tmp*(v_na-v_tmp)+ ...
        g_k*n_tmp^4*(v_k-v_tmp)+g_l*(v_l-v_tmp)+i_ext_i)/c;
    m_inc=alpha_m(v_tmp)*(1-m_tmp)-beta_m(v_tmp)*m_tmp;
    h_inc=alpha_h(v_tmp)*(1-h_tmp)-beta_h(v_tmp)*h_tmp;
    n_inc=alpha_n(v_tmp)*(1-n_tmp)-beta_n(v_tmp)*n_tmp;
    
    v_i(k+1)=v_i(k)+dt*v_inc - w(k) ;
    m(k+1)=m(k)+dt*m_inc;
    h(k+1)=h(k)+dt*h_inc;
    n(k+1)=n(k)+dt*n_inc;

    if v_i(k) >= -50
        w(k+1)= w(k) + dt*(a_k*(v_i(k) - v_l) - w(k) + b_k*t_k)/t_k;
        if v_i(k+1) < -50
            i_counter = i_counter+1;
            spiketimes_i(i_counter) = k;
        end

    else
        w(k+1)= w(k) + dt*(a_k*(v_i(k) - v_l) - w(k))/t_k;
    end
    
end

t_i=(0:m_steps)*dt;

w_i = w;

end