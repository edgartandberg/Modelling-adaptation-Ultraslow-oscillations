function [u,a,spikecount] = Synaptic_input(previous_chain,u,t_ref,u_rest,num_end,u_th,u_spike,num_start,u_hp,a_gain,weight)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



a = t_ref;
spikecount = zeros(1, num_end - 2 - a);


u = previous_chain;


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