function [u,a,spikecount] = LIF(u_rest, R, Iext, t_m, t_ref, num_end, u_th, u_spike, u_hp)

% LiF model as ODE
% ------------------------
% This function calculates voltage from solved diff. equation for the leaky
% integrate and fire model (du/dt).
% Loops over time steps, spikes if threshold is reached, before
% hyperpolarizing and returning to resting potential
% Also makes 'spikecount' variable for raster plot, by adding a 1 to the
% empty array wehenever a spike occurs
% ------------------------


dt=0.25; %in seconds
spk_times=[];
counter=0;
u(2)=u_rest; % intitial conditions u

%du/dt = (-u + u_rest + R*Iext)/tau_m

t=2;

while t<=800
    u(t)= u(t-1) + dt*(-(u(t-1) - u_rest) + R.*Iext(t-1))/t_m;
    if (u(t)>u_th)
        u(t) = u_spike;
        u(t+1) = u_hp;
        u(t+2) = u_rest;
        t = t+2;
        counter=counter+1;
        spk_times(counter)=t;
    end
    t = t+1;
end


a = t_ref;
spikecount = zeros(1, num_end - 2 - a);

