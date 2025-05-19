function [RI, Iext, noise] = Gen_Current(R, Iext,num_start, num_end, num_pulses, I_0)
% Takes in input current w parameters and spits out function for
% current
% ------------------------
%   takes in resistance, single value for input current, parameters for how
%   long and how many pulses. Gives us function for current that is used in
%   solution of diff. equation for u, with added noise
% ------------------------

noise = 8*1e-8*randn(size(Iext)); % noise when current is injected
num_width = round(num_end - num_start);
 


freq = 0.015 * num_pulses; % adjusting frequency for right amount of pulses

time = num_start:1:num_end;
signal =  I_0*sin(time*freq); % modelling input current as sine function for oscillations


Iext(num_start : num_start + num_width) = I_0; 
% change to Iext = I_0 for step input, Iext = signal for osc input

%Iext = Iext + noise; %add this to add noise

Iext(1:num_start) = 0; % resetting signal before and after current to 0
Iext(num_end:end) = 0;


RI = R.*Iext;
end