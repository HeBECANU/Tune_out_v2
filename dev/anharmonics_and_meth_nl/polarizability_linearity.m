fsamp=linspace(1,1000,1e5)*1e12;
polz_au=aprox_he_polz(fsamp);
stfig('atomic polz');
clf
plot(fsamp*1e-12,polz_au)
ylim([-200,600])
yline(0)

% find tune-out using crossing search
fminopt = optimset('TolX',1e-20,'TolFun',1e-6);   %,'Display','iter'
tune_out_freq=fzero(@(x) aprox_he_polz(x),f2wl(413e-9),fminopt);  
au_polz_at_to=aprox_he_polz(tune_out_freq);
f2wl(tune_out_freq)
hold on
plot(tune_out_freq*1e-12,au_polz_at_to,'xr')
hold off
fprintf('found tune out at %.6f THz, %.6f nm, polz %g au \n',tune_out_freq*1e-12,f2wl(tune_out_freq)*1e9,aprox_he_polz(tune_out_freq))

% find tune-out using minimization of abs value
% fminopt = optimset('TolX',1,'TolFun',1e-12);  
% tune_out_freq=fminsearch(@(x) abs(aprox_he_polz(x)),600e12,fminopt);  
% hold on
% plot(tune_out_freq*1e-12,aprox_he_polz(tune_out_freq),'xr')
% hold off
% fprintf('found tune out at %.6f THz, %.6f nm, polz %g au \n',tune_out_freq*1e-12,f2wl(tune_out_freq)*1e9,aprox_he_polz(tune_out_freq))
% 

%% derivatives in AU/hz
to_polz_au_deriv=arrayfun(@(y) derivest(@(x) aprox_he_polz(x),tune_out_freq,'DerivativeOrder',y),1:4);
fprintf('au polz deriv /hz')
to_polz_au_deriv(1)
to_polz_au_deriv(2)
%%
to_polz_si_deriv=arrayfun(@(y) derivest(@(x) outselect(2,@aprox_he_polz,x),tune_out_freq,'DerivativeOrder',y,'FixedStep',0.1),1:4);
to_polz_si_deriv(1)
to_polz_si_deriv(2)

%% try fiting a linear about the tune out

fsamp=col_vec(linspace(-1,1,1e5))*4e10+tune_out_freq;
[polz_au,~]=aprox_he_polz(fsamp);
stfig('atomic polz');
clf
xscale=1e-9;
%yscale=1e44;
yscale=1;
plot((fsamp-tune_out_freq)*xscale,polz_au*yscale)
xlabel('$\omega_{\mathrm{TO}}-\omega$ ($2\pi$ GHz)')
ylabel(sprintf('$\\alpha$ ($10^{%g} \\mathrm{C}\\cdot \\mathrm{m}^{2} \\cdot\\mathrm{V}^{-1}$)',-log10(yscale)))
yline(0)


predictor=(fsamp-tune_out_freq)*1e-9;
response=polz_au*yscale;

fit_in=cat(2,predictor,response,response*nan);
meth_lin_fit=fit_poly_with_int(fit_in,1,0,0);
fprintf('lin fit intercept %f MHz \n',(meth_lin_fit.x_intercept.val/xscale)*1e-6)


fit_in=cat(2,predictor,response,response*nan);
meth_lin_fit=fit_poly_with_int(fit_in,2,0,0);
fprintf('quad fit intercept %f MHz \n',(meth_lin_fit.x_intercept.val/xscale)*1e-6)