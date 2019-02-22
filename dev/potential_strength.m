hebec_constants

e0 = const.epsilon0;
c = const.c;
a0 = const.a0;
TO_m = 413.085*1e-9;
TO_freq = c/TO_nm;
unc_freq = 700/sqrt(1000)*1e6;%Hz
unc_m = ((TO_m)^2/c)*unc_freq;
%unc_m = abs(c/(TO_freq+unc_freq)-TO_nm);
alpha = unc_m*1.913*a0^3*1e9;%in cgs m^3
alpha = alpha*4*pi*e0;%in SI Fm^2
P_beam = 100*1e-3; %Watts, average power of beam
w_beam = 10*1e-6;%m
I_beam = 2*P_beam/(pi*w_beam^2);

U = (1/2)*(1/(e0*c))*alpha*I_beam
