%folowing http://iopscience.iop.org/article/10.1209/epl/i2006-10284-4

hebec_constants
%rabi freq
a=9.1e-8;
lambda=1557e-9;
transtion_freq=2*pi*const.c/lambda;
power=2;
radius=280e-6;
omega=sqrt((2*pi*(const.c^2)/(const.h*(transtion_freq^3)))*a*power/(radius^2));
fprintf('transtion lifetime 1/rate %.2e \n',1/omega)

%%
%rabi freq
a=6e-9;
lambda=427e-9;
transtion_freq=2*pi*const.c/lambda;
power=50e-3;
radius=50e-6;
omega=sqrt((2*pi*(const.c^2)/(const.h*(transtion_freq^3)))*a*power/(radius^2))
fprintf('transtion lifetime 1/rate %.2e \n',1/omega)


%% calculating wavelength
