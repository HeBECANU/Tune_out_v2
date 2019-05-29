%folowing http://iopscience.iop.org/article/10.1209/epl/i2006-10284-4

hebec_constants
%rabi freq
A=9.1e-8;
lambda=1557e-9;
transition_freq=2*pi*const.c/lambda;
power=2;
radius=136e-6;
D= 10e-3;
Intensity = 2*power/radius^2;
% omega=2*sqrt((2*pi*(const.c^2)/(const.hb*(transtion_freq^3)))*a*power/(radius^2));
omega=sqrt((2*pi*const.c^2)/(const.hb*transition_freq^3)*A*Intensity);
disp('-- 1557nm --')
fprintf('transtion rate            %.2e \n',omega)
fprintf('transtion lifetime 1/rate %.2e \n',1/omega)

%%
%rabi freq for 427
A=6.4e-9;
lambda=427e-9;
transtion_freq=2*pi*const.c/lambda;
power=50e-3;
radius=20e-6;
intensity = 2*power/radius^2;
omega_427=sqrt((2*pi*const.c^2)/(const.hb*transition_freq^3)*A*intensity);
disp('-- 427nm --')
fprintf('transtion rate            %.2e \n',omega_427)
fprintf('transtion lifetime 1/rate %.2e \n',1/omega_427)


%%
% Dumb comparison to 51D2
A=8e2; %Using the value ffor the 3P1...
lambda=402.74-9;
transtion_freq=2*pi*const.c/lambda;
power=20e-3;
radius=1e-3;
intensity = 2*power/radius^2;
omega_427=sqrt((2*pi*const.c^2)/(const.hb*transition_freq^3)*A*intensity);
disp('-- 402nm --')
fprintf('transtion rate            %.2e \n',omega_427)
fprintf('transtion lifetime 1/rate %.2e \n',1/omega_427)