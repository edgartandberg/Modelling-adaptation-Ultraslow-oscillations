clc

dydt = @(t,y) 2*t;
ode45(dydt, [-1,1],1)

%%

N = 1; % number of neurons

tau_m = 10;  % membrane time constant
v_r = 0;  % reset potential
v_th = 15;  % threshold potential
I_c = 20; % input current
R = 10^7; % resistance


t = 0:1:1000;

v =  R*I_c*(1-exp(-t/tau_m));

plot(t,v)

%%

dvdt = @(t,v) (R*I_c-v)/tau_m;

ode45(dvdt, [0,1000],0)
