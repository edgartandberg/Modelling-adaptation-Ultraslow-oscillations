function alpha_n=alpha_n(v);
alpha_n=0.01*(-45.0d0-v)./(exp((-45-v)/10)-1);
