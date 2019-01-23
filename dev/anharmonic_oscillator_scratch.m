<<<<<<< HEAD
% clear all

%% Parameters
X0 = [1 0];
trap_freq = 155*2*pi;
T_win = linspace(0,9,10000);
P = [1 trap_freq^2 0 0 0];%[m, spring, damping, cubic, quartic] 
% P=P0.*(1+0.05*randn(1,5)); 

t_offset = min(T_win):5:max(T_win);
t_width = 50;

beta0 = [sqrt(P(2)/P(1))/(2*pi), 0,0,0];
f_samp = 9.8;
T_samp = 1/f_samp;
N_samp = 150;
N_zone = 2*floor(trap_freq/f_samp);
T_det = linspace(0,N_samp*T_samp,N_samp);
%N_zone = trap_freq/(f_samp/2);
f_guess = 2*mod(trap_freq,f_samp);
P_downsamp = [f_guess/(2*pi) 0 0 0]; %amp freq damping offset
fit_coef_names = {'Freq', 'phase', 'damp', 'offset'};



%% Heavy lifting

[t,y] = ode45(@(t,X) trap_DE(t,X,P),T_win,X0);
%[t,y_apx] = ode45(@(t,X) trap_DE(t,X,P0),T_win,X0);

data = [t,y(:,1)];


fit = fitnlm(t,y(:,1),@damped_sine,beta0,'CoefficientNames',fit_coef_names)
proj = down_sample(@damped_sine,fit.Coefficients.Estimate,T_det);
fit_downsamp = fitnlm(T_det,proj,@damped_sine, P_downsamp,'CoefficientNames',fit_coef_names)
T_fit = linspace(0,max(T_det),500);
5*f_samp+abs(fit_downsamp.Coefficients.Estimate(1))
=======
%demonstration of imprefect fit when the trap is not prefectly harmonic
% trying to come up with some way to add extra terms to the fit to account for it
% at the 
% clear all
%

%% Parameters
inital_pos_vel = [1 0]; %inital position,velocity
trap_freq = 500*2*pi;
damping=0;%1e-0;
t_samp_solver = linspace(0,2,1e6);
mass=1;
spring_const=mass*trap_freq^2;
damping_ratio=damping/(2*sqrt(mass*spring_const));
damping_time=1/(trap_freq*damping_ratio);
resonant_freq=trap_freq*sqrt(1-2*damping_ratio^2);
cubic_trap=1e5;
quartic_trap=0;
P = [1 spring_const damping cubic_trap quartic_trap];%[m, spring, damping, cubic, quartic] 
% P=P0.*(1+0.05*randn(1,5)); 

t_offset = min(t_samp_solver):5:max(t_samp_solver);
t_width = 50;

fit_guess_fundemental = [resonant_freq, 0,P(3),0]; %'Freq', 'phase', 'damp', 'offset'
fit_guess_2nd= [resonant_freq, 0,P(3),0,1e-3,0];

P_downsamp = [f_guess/(2*pi) 0 0 0]; %amp freq damping offset




% Heavy lifting

[t_ode_out,y_ode_out] = ode45(@(t,X) trap_DE(t,X,P),t_samp_solver,inital_pos_vel);
%[t,y_apx] = ode45(@(t,X) trap_DE(t,X,P0),T_win,X0);

%%
data = [t_ode_out,y_ode_out(:,1)];
opts = statset('TolFun',1e-5,'TolX',1e-5);
fit_coef_names = {'Freq', 'phase', 'damp', 'offset'};
fit_fundemental = fitnlm(t_ode_out,y_ode_out(:,1),@damped_sine_fundemental,fit_guess,'CoefficientNames',fit_coef_names,'Options',opts);
%%
fit_coef_names = {'Freq', 'phase', 'damp', 'offset','harmonic_amp','harmonic_phase'};
fit_2nd = fitnlm(t_ode_out,y_ode_out(:,1),@damped_sine_2nd,fit_guess_2nd,'CoefficientNames',fit_coef_names,'Options',opts);

%%
% Graphical output
figure(2)
clf
subplot(3,1,1)
plot(t_ode_out,y_ode_out(:,1),'k-')
hold on
plot(t_ode_out,1e-3+predict(fit_fundemental,t_ode_out),'r-');
plot(t_ode_out,1e-3+predict(fit_2nd,t_ode_out),'b-');
title('raw motion')
hold off

%residuals
subplot(3,1,2)
fit_residuals_fund=y_ode_out(:,1)-predict(fit_fundemental,t_ode_out);
fit_residuals_2nd=y_ode_out(:,1)-predict(fit_2nd,t_ode_out);
plot(t_ode_out,fit_residuals_fund,'k-')
hold on
plot(t_ode_out,fit_residuals_2nd,'r-')
hold off

subplot(3,1,3)
fft_resid_fund=fft_tx(t_ode_out,fit_residuals_fund,'padding',10,'window','hamming');
fft_resid_2nd=fft_tx(t_ode_out,fit_residuals_2nd,'padding',10,'window','hamming');
plot(fft_resid_fund(1,:),abs(fft_resid_fund(2,:)),'k-')
hold on
plot(fft_resid_2nd(1,:),abs(fft_resid_2nd(2,:)),'r-')
hold off
xlim([0,3*trap_freq/2*pi])




%%


%proj = down_sample(@damped_sine,fit.Coefficients.Estimate,T_det);
%fit_downsamp = fitnlm(T_det,proj,@damped_sine, P_downsamp,'CoefficientNames',fit_coef_names)
%T_fit = linspace(0,max(T_det),500);
%5*f_samp+abs(fit_downsamp.Coefficients.Estimate(1))
>>>>>>> master

%N_zone-mod(N_zone,2)*f_samp/2 -(-1)^mod(floor(N_zone),2)*fit_downsamp.Coefficients.Estimate(1)

% 
% Y = proj;
% delay_freq = zeros(N_samp,1);
% for i=1:numel(t_offset)-1
%     t0 = t_offset(i);
%     t_sel = t(t>t0 & t<t0+t_width)-t_offset(i);
%     amp_sel = Y(t>t0 & t<t0+t_width);
%     fit_sel = fitnlm(t_sel,amp_sel,@damped_sine,beta0);
%     delay_freq(i) = fit_sel.Coefficients.Estimate(1);
% end

<<<<<<< HEAD
%% Graphical output

subplot(3,1,1)
plot(t,y(:,1),'b')
hold on
plot(t,damped_sine(fit.Coefficients.Estimate,t),'r.');
title('raw motion')
hold off

subplot(3,1,2)
plot(T_det,proj,'x')
hold on
plot(T_fit,damped_sine(fit_downsamp.Coefficients.Estimate,T_fit))
hold off
title('Downsampled signal')
subplot(3,1,3)




function y = damped_sine(p,t)
    y = exp(-abs(p(3))*t).*sin(2*pi*p(1)*t-p(2))+p(4);
end

=======
% f_samp = 9.8;
% T_samp = 1/f_samp;
% N_samp = 150;
% N_zone = 2*floor(trap_freq/f_samp);
% T_det = linspace(0,N_samp*T_samp,N_samp);
% %N_zone = trap_freq/(f_samp/2);
% f_guess = 2*mod(trap_freq,f_samp);


% subplot(3,1,2)
% plot(T_det,proj,'x')
% hold on
% plot(T_fit,damped_sine(fit_downsamp.Coefficients.Estimate,T_fit))
% hold off
% title('Downsampled signal')


function y = damped_sine_fundemental(p,t)
    y = exp(-abs(p(3))*t).*sin(2*pi*p(1)*t-p(2))+p(4);
end

function y = damped_sine_2nd(p,t)
    amp=exp(-abs(p(3))*t);
    y = amp.*sin(2*pi*p(1)*t-p(2))+p(4)+p(5)*amp.*sin(4*pi*p(1)*t-p(6));
end

>>>>>>> master
% function y = mod_sine(p,t)
%     y = exp(-p(3)*t).*sin(p(1).*t-p(2))+p(4);
% end

function dX = trap_DE(t,X,p)
<<<<<<< HEAD
    dX(1) = X(2);% - p(3)*X(1);
=======
    dX(1) = X(2); %- p(3)*X(1); %comment out after the first term to turn off damping
>>>>>>> master
    dX(2) = (-p(2)*X(1) - (p(4)*X(1)^2)/2 - (p(5)*X(1)^3)/6)/p(1);
    dX = dX';
end %function

function f = objfun(P,data,IC)
    [t,y_apx] = ode45(@(t,X) trap_DE(t,X,P),data(:,1),IC);
    f = sum((data(:,2)-y_apx(:,1)).^2); %L2 error
end

function P = down_sample(func, P,T)
    P = func(P,T);
end

