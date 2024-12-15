tic

% Parameters
t_k = 400; % time constant for adaptation current, indexing over k neurons
a_k = 0.5; % adaptation coupling
b_k = 5000; % adaptation gain
u_rest = -70; % resting potential
t_m = 80; % time point for spike
u_r = -55; % voltage reset
u_th= 50;
dt=0.01; %in seconds
spk_times=[];
counter=0;
u(1)=40;
w(1)=3;
I(1,1:25000)=0;
I(1,25000:75000)=200;
I(1,75000:100000)=0;

for t=2:100000
    u(t)= u(t-1) + dt*(-(u(t-1) - u_rest) - 1*w(t-1) + I(t-1))/t_m;
    if (u(t)>=u_th)
        u(t)=u_r;
        counter=counter+1;
        spk_times(counter+1)=t;
        w(t)= w(t-1) + dt*(a_k*(u(t-1) - u_rest) -w(t-1) + b_k*t_k)/t_k;
    else
        %w(t)= w(t-1) + dt*(a_k*(u(t-1) - u_rest) -w(t-1) + b_k*t_k*length(spk_times))/t_k;
        w(t)= w(t-1) + dt*(a_k*(u(t-1) - u_rest) -w(t-1))/t_k;
    end
end
figure
subplot(3,1,1);
plot(u,'-*')
subplot(3,1,2);
plot(I)
subplot(3,1,3);
plot(w)

toc