function [v_i, t_i, i_counter, m, h, n, spiketimes_i, i_ext_i_all, z, b] = I_Sompolinsky_adaptation(v_i, k, t_final,dt,i_ext_i, m_i, h_i, n_i, z_i, b_i)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

dt05=dt/2;


c = 1; % Capacitance
g_na = 120; % Conductance for sodium
g_k = 20; % Conductance for potassium
g_a = 40;  % A-current conductance
g_nap = 0.2; % Persistent sodium current conductance 
g_l = 0.1; % Leak conductance
g_z = 0;  % slow potassium conductance / causes adaptation
v_na = 55; % Sodium reversal potential
v_k = -70; % Potassium reversal potential
v_a = -75; % A-current reversal potential
v_l = -65; % Leak reversal potential
tau_z = 60; % slow potassium current time constant
tau_b = 10; % A-current time constant

m=m_i;
h=h_i;
n=n_i;
z=z_i;
b=b_i;

spiketimes_i = [];
m_steps=round(t_final/dt);

if k < 2
i_ext_i = i_ext_i(k-1);
end


v_i(1)=-70; % potential for excitatory neuron, initial condition

% gating variables, dimensionless
m(1)=alpha_m(v_i(1))/(alpha_m(v_i(1))+beta_m(v_i(1)));
h(1)=0.5; 
n(1)=0.5;
z(1)=0.0;
a(1)=0.0;
b(1)=0.0;
s(1)=0.0;


i_counter=0;
i_ext_i_all = zeros(1,m_steps+1);




    %i_ext_i = 1.5; % in micro A

    v_inc=(g_na*m(k-1)^3*h(k-1)*(v_na-v_i(k-1))+ g_k*n(k-1)^4*(v_k-v_i(k-1))+g_l*(v_l-v_i(k-1)) ...
        + g_a*a_infty(v_i(k-1))*b(k-1)*(v_a-v_i(k-1)) ...
        + g_nap*s_infty(v_i(k-1))*(v_na-v_i(k-1)) + i_ext_i)/c;

    m_inc=alpha_m(v_i(k-1))*(1-m(k-1))-beta_m(v_i(k-1))*m(k-1);
    h_inc=alpha_h(v_i(k-1))*(1-h(k-1))-beta_h(v_i(k-1))*h(k-1);
    n_inc=alpha_n(v_i(k-1))*(1-n(k-1))-beta_n(v_i(k-1))*n(k-1);
    %z_inc=(z_infty(v_i(k-1))-z(k-1))/tau_z;
    b_inc=(b_infty(v_i(k-1))-b(k-1))/tau_b;


    v_tmp=v_i(k-1)+dt05*v_inc;
    m_tmp=m(k-1)+dt05*m_inc;
    h_tmp=h(k-1)+dt05*h_inc;
    n_tmp=n(k-1)+dt05*n_inc;
    %z_tmp=z(k-1)+dt05*z_inc;
    b_tmp=b(k-1)+dt05*b_inc;


    v_inc=(g_na*m_tmp^3*h_tmp*(v_na-v_tmp)+ ...
        g_k*n_tmp^4*(v_k-v_tmp)+g_l*(v_l-v_tmp) + g_a*a_infty(v_tmp)*b(k-1).*(v_a-v_tmp) ...
        + g_nap*s_infty(v_tmp).*(v_na-v_tmp) + i_ext_i)/c;


    m_inc=alpha_m(v_tmp)*(1-m_tmp)-beta_m(v_tmp)*m_tmp;
    h_inc=alpha_h(v_tmp)*(1-h_tmp)-beta_h(v_tmp)*h_tmp;
    n_inc=alpha_n(v_tmp)*(1-n_tmp)-beta_n(v_tmp)*n_tmp;
    %z_inc=((z_infty(v_tmp)-z(k-1))/tau_z)*z_tmp;
    b_inc=(b_infty(v_tmp)-b(k-1))/tau_b;


    if length(v_inc) > 1

    v_i(k)=v_i(k-1)+dt*v_inc(k-1);
    m(k)=m(k-1)+dt*m_inc(k-1);
    h(k)=h(k-1)+dt*h_inc(k-1);
    n(k)=n(k-1)+dt*n_inc(k-1);
    %z(k)=z(k-1)+dt*z_inc(k-1);
    b(k)=b(k-1)+dt*b_inc(k-1);

    else

    v_i(k)=v_i(k-1)+dt*v_inc;
    m(k)=m(k-1)+dt*m_inc;
    h(k)=h(k-1)+dt*h_inc;
    n(k)=n(k-1)+dt*n_inc;
    %z(k)=z(k-1)+dt*z_inc;
    b(k)=b(k-1)+dt*b_inc;

    end


    if v_i(k-1) > 20
        if v_i(k) < 20
            i_counter = i_counter+1;
            spiketimes_i(i_counter) = k;
        end
    end

    %i_ext_i_all(k-1) = i_ext_i;    

t_i=(0:m_steps)*dt;


end