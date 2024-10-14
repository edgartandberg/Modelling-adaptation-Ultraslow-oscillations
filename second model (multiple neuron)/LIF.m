function [u,a,spikecount] = LIF(u_rest, RI, t, tau_m, t_ref, num_end, u_th, u_spike, u_hp, a_gain, num_start, connectivitymatrix)

% LiF model as ODE
% ------------------------
% This function calculates voltage from solved diff. equation for the leaky
% integrate and fire model (du/dt).
% Loops over time steps, spikes if threshold is reached, before
% hyperpolarizing and returning to resting potential
% Also makes 'spikecount' variable for raster plot, by adding a 1 to the
% empty array wehenever a spike occurs
% ------------------------

u = u_rest + RI.*(1-exp(-(t/tau_m))); %solution for above equation for u at u(0) = u_rest



a = t_ref;
spikecount = zeros(1, num_end - 2 - a);

for ii = 1 : num_end - 2 - a % Adjust loop range to avoid out-of-bounds errors
    if u(ii:ii+1000) < u_th % if below threshold for 1000 ms, resets adaptation
        a = t_ref;
    end

    if u(ii + 1) > u_th
        u(ii) = u_spike;
        spikecount(ii-num_start) = 1;
        u(ii + 1) = u_hp;
        u(ii + 2 : ii + 2 + a) = u_rest;
        a = round(a + a_gain); % Increment by adaptation gain
    end
end


end


