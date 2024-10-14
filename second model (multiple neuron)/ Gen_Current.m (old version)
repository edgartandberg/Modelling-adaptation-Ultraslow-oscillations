function [RI, Iext] = Gen_Current(R, Iext,num_start, num_end, num_pulses, I_0)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here


noise = 5*1e-8*randn(size(Iext)); % noise when current is injected
num_width = round(num_end - num_start);
 


freq = 0.001 * num_pulses; % adjusting frequency for right amount of pulses

time = num_start:1:num_end;
signal =  I_0*sin(time*freq);


Iext(num_start : num_start + num_width) = signal; % modelling input current as sine function for oscillations
% change to = I_0 for step input

Iext = Iext + noise;
Iext(1:num_start) = 0; % resetting signal before and after current to 0
Iext(num_end:end) = 0;


RI = R.*Iext;
end
