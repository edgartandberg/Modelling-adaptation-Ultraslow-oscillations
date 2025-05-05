%% slow potassium current


clc
clearvars

t_final=2000; % in ms
dt=0.01;
dt05=dt/2;


c = 1; % Capacitance
g_na = 120; % Conductance for sodium
g_k = 10; % Conductance for potassium
g_a = 60;  % A-current conductance
g_nap = 0.5; % Persistent sodium current conductance 
g_l = 0.1; % Leak conductance
g_z = 10;  % slow potassium conductance / causes adaptation
v_na = 55; % Sodium reversal potential
v_k = -70; % Potassium reversal potential
v_a = -75; % A-current reversal potential
v_l = -65; % Leak reversal potential
tau_z = 60; % slow potassium current time constant
tau_b = 10; % A-current time constant




m_steps=round(t_final/dt);



v_e(1)=-70; % potential for excitatory neuron, initial condition

% gating variables, dimensionless
m(1)=alpha_m(v_e(1))/(alpha_m(v_e(1))+beta_m(v_e(1)));
h(1)=0.5; 
n(1)=0.5;
z(1)=0.0;
a(1)=0.0;
b(1)=0.0;
s(1)=0.0;


e_counter=0;
i_ext_e_all = zeros(1,m_steps+1);



for k=1:m_steps

    i_ext_e = 1.5; % in micro A

    v_inc=(g_na*m(k)^3*h(k)*(v_na-v_e(k))+ g_k*n(k)^4*(v_k-v_e(k))+g_l*(v_l-v_e(k)) ...
        + g_z*z(k)*(v_k-v_e(k)) + g_a*a_infty(v_e(k))*b(k)*(v_a-v_e(k)) ...
        + g_nap*s_infty(v_e(k))*(v_na-v_e(k)) + i_ext_e)/c;

    m_inc=alpha_m(v_e(k))*(1-m(k))-beta_m(v_e(k))*m(k);
    h_inc=alpha_h(v_e(k))*(1-h(k))-beta_h(v_e(k))*h(k);
    n_inc=alpha_n(v_e(k))*(1-n(k))-beta_n(v_e(k))*n(k);
    z_inc=(z_infty(v_e(k))-z(k))/tau_z;
    b_inc=(b_infty(v_e(k))-b(k))/tau_b;


    v_tmp=v_e(k)+dt05*v_inc;
    m_tmp=m(k)+dt05*m_inc;
    h_tmp=h(k)+dt05*h_inc;
    n_tmp=n(k)+dt05*n_inc;
    z_tmp=z(k)+dt05*z_inc;
    b_tmp=b(k)+dt05*b_inc;


    v_inc=(g_na*m_tmp^3*h_tmp*(v_na-v_tmp)+ ...
        g_k*n_tmp^4*(v_k-v_tmp)+g_l*(v_l-v_tmp) + g_z*z(k)*(v_k-v_tmp) + g_a*a_infty(v_tmp)*b(k)*(v_a-v_tmp) ...
        + g_nap*s_infty(v_tmp)*(v_na-v_tmp) + i_ext_e)/c;


    m_inc=alpha_m(v_tmp)*(1-m_tmp)-beta_m(v_tmp)*m_tmp;
    h_inc=alpha_h(v_tmp)*(1-h_tmp)-beta_h(v_tmp)*h_tmp;
    n_inc=alpha_n(v_tmp)*(1-n_tmp)-beta_n(v_tmp)*n_tmp;
    z_inc=((z_infty(v_tmp)-z(k))/tau_z)*z_tmp;
    b_inc=(b_infty(v_tmp)-b(k))/tau_b;


    v_e(k+1)=v_e(k)+dt*v_inc;
    m(k+1)=m(k)+dt*m_inc;
    h(k+1)=h(k)+dt*h_inc;
    n(k+1)=n(k)+dt*n_inc;
    z(k+1)=z(k)+dt*z_inc;
    b(k+1)=b(k)+dt*b_inc;



    if v_e(k) > 20
        if v_e(k+1) < 20
            e_counter = e_counter+1;
            spiketimes_e(e_counter) = k;
        end
    end

    i_ext_e_all(k) = i_ext_e;    
end




%plot(v_e)
plot(diff(spiketimes_e))
