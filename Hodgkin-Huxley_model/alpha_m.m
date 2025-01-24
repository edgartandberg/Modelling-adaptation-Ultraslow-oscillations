function alpha_m=alpha_m(v);
if abs(v+30)>1.0e-8,
    alpha_m=(v+30)/10./(1-exp(-(v+30)/10));
else
    alpha_m=1;
end;