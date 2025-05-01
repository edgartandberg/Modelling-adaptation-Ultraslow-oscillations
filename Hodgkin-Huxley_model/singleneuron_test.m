clc
clearvars

tic

t_final=100; % in ms
dt=0.01;
dt05=dt/2;
m_steps=round(t_final/dt);


i_ext_e = 3; % in micro A
i_ext_i = 3; % in micro A
v_e = zeros(1, m_steps);
v_i = zeros(1, m_steps);
v_e(1)=-70; % initial condition
v_i(1)=-70; % initial condition
size_network = 1;

% % initialize arrays for gating variables, excitatory
m_e = zeros(1,m_steps);
h_e = zeros(1,m_steps);
n_e = zeros(1,m_steps);

% inhibitory
m_i = zeros(1,m_steps);
h_i = zeros(1,m_steps);
n_i = zeros(1,m_steps);

% intitialize arrays for conductances, excitatory
g_e_all = cell(1, size_network);
I_syn_e_all = zeros(size_network, m_steps);

% inhibitory
g_i_all = cell(1, size_network);
I_syn_i_all = zeros(size_network, m_steps);


% initital conditions for gating variables, excitatory
m_e(1)=alpha_m(v_e(1))/(alpha_m(v_e(1))+beta_m(v_e(1)));
h_e(1)=0.5; 
n_e(1)=0.5;

% inhibitory
m_i(1)=alpha_m(v_i(1))/(alpha_m(v_i(1))+beta_m(v_i(1)));
h_i(1)=0.5; 
n_i(1)=0.5; 

g_e = zeros(1, m_steps); % intialize synaptic conductance
g_i = zeros(1, m_steps);

for k=2:m_steps
    for  j = 1:size_network



    [v_e, t_e, e_counter, m_e, h_e, n_e, spk_times_e] = excitatory_HH_pass2(v_e, k, t_final,dt,i_ext_e, m_e, h_e, n_e);
    [v_i, t_i, i_counter, m_i, h_i, n_i, spk_times_i] = inhibitory_HH_pass2(v_i, k, t_final,dt,i_ext_i, m_i, h_i, n_i);
    

    % spk_array_e(k) = spk_times_e;
    % spk_array_i(k) = spk_times_i;

    g_e_old = g_e(k-1);
    g_i_old = g_i(k-1);

   

    % I_syn_e_all(j, k) = I_syn_e; % Store the synaptic current for each network
    % I_syn_i_all(j, k) = I_syn_i; % Store the inhibitory synaptic current 
    spiketrains_e_all{j} = v_e;
    spiketrains_i_all{j} = v_i;

    i_ext_e_all(j,k) = i_ext_e;
    i_ext_i_all(j,k) = i_ext_i;


    % output the spike trains from the excitatory neurons
    m_e_all{j} = m_e;
    h_e_all{j} = h_e;
    n_e_all{j} = n_e;

    m_i_all{j} = m_i;
    h_i_all{j} = h_i;
    n_i_all{j} = n_i;

    end


end
