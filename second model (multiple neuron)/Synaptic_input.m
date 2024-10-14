function [u,a,spikecount] = Synaptic_input(previous_chain,u,t_ref,u_rest,num_end,u_th,u_spike,current_start,u_hp,a_gain,weight)
% Synaptic input for generating sequence
% -------------------------
% Takes in previous chain from u_all with weight from conn.matrix 
% for generating the next chain in sequence
% -------------------------



a = t_ref;
spikecount = zeros(1, num_end - 2 - a);


u = previous_chain*weight;
u = [zeros(1, current_start*0.1), u(1:end-current_start*0.1)];


for ii = current_start: num_end - 2 - a % Adjust loop range to avoid out-of-bounds errors
    if u(ii:ii+1000) < u_th % if below threshold for 1000 ms, resets adaptation
        a = t_ref;
    end

    if u(ii + 1) > u_th
        u(ii) = u_spike;
        spikecount(ii-current_start) = 1;
        u(ii + 1) = u_hp;
        u(ii + 2 : ii + 2 + a) = u_rest;
        a = round(a + a_gain); % Increment by adaptation gain
    end
end

end