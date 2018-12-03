%User config
num_shot = 30;
pulse_range = 6:2:40;
%Note RF power drops after 20 pulses, so we'll only fit these
pow_fit_mdl = @(p,x) p(1) * p(2).^(x*p(3));
coefs = [1.2e3,0.8,0.2];
trap_freqs = data.osc_fit.model_coefs(:,2,1);


%Setting up
all_N = data.mcp_tdc.num_counts; % idx shot #
pulse_N = data.mcp_tdc.al_pulses.num_counts; %idx (shot #, pulse #)
% kill nans
pulse_mask = ~isnan(sum(pulse_N,2));
freq_mask = ~isnan(trap_freqs);
mask = pulse_mask & freq_mask;
trap_freqs = trap_freqs(mask);
all_N = all_N(mask);
pulse_N = pulse_N(mask,:);

figure()
plot(normalize(all_N),normalize(trap_freqs),'x');
xlabel('Normalized number')
ylabel('Normalized fit freq')


%%
%Trim down to user config sizes
pulse_idx = 1:size(pulse_N,2);
all_N = all_N(1:num_shot);
trap_freqs = trap_freqs(1:num_shot);
fit_idx = pulse_idx(1:max(pulse_range));
fit_num = pulse_N(:,1:max(pulse_range));



close all
figure()

num_pulses = numel(pulse_range);
scaling_data = cell(num_pulses);
coef_all = zeros(num_pulses,num_shot,3);
frac_errs = zeros(num_pulses,num_shot);
rem_pops = zeros(num_pulses,num_shot);
rem_errs = zeros(num_pulses,num_shot);
dropped = zeros(num_pulses,num_shot);
disp('Fitting!')
for m=1:num_pulses
    fit_idx_m = fit_idx(1:pulse_range(m));
    for n=1:num_shot
        fit_dat_n = fit_num(n,1:pulse_range(m));
        pow_fit = fitnlm(fit_idx_m,fit_dat_n,pow_fit_mdl,coefs);
        coef_all(m,n,:) = pow_fit.Coefficients.Estimate;
        dropped(m,n) = sum(pulse_N(n,1:pulse_range(m)));
        rem_pops(m,n) = all_N(n) - sum(pulse_N(n,1:pulse_range(m)));
        rem_est = (coef_all(m,n,1)./(1-coef_all(m,n,2).^coef_all(m,n,3)))-sum(pulse_N(n,1:pulse_range(m)));
        rem_errs(m,n) = (rem_pops(m,n) - rem_est)/rem_pops(m,n);
    end   
    subplot(2,3,1)
    histogram(coef_all(m,:,2),linspace(0.8,0.9,20));
    hold on
    subplot(2,3,2)
    histogram(coef_all(m,:,3),linspace(0.1,0.2,20));
    hold on
    subplot(2,3,3)
    recon_N = coef_all(m,:,1)./(1-coef_all(m,:,2).^coef_all(m,:,3));
    frac_err = (all_N-recon_N)./all_N;
    histogram(frac_err,linspace(-0.4,0.4,20))
    hold on
    frac_errs(m,:) = frac_err;
    
    
end
subplot(2,3,1)
title('Decay coefficient')
subplot(2,3,2)
title('Decay power')
subplot(2,3,3)
title('Fractional error in N')

subplot(2,3,4)
errorbar(pulse_range, mean(frac_errs,2),std(frac_errs'))
xlim([0,max(pulse_range)+1])
ylim([-1,1])
xlabel('Pulse number')
ylabel('Error in N')
title('Error scaling')
subplot(2,3,5)
X = mean(dropped./all_N,2);
Y=abs(mean(frac_errs,2));
Y_err = std(frac_errs,[],2);
errorbar(X,log(Y),0.434*Y_err./Y,'x')
% set(gca,'Yscale','log')
xlabel('Outcoupled fraction')
ylabel('log (error)')
title('lost fraction vs precision')
subplot(2,3,6)
X = mean(dropped./all_N,2);
X_err = std(dropped./all_N,[],2);
Y = mean(rem_errs,2);
Y_err = std(rem_errs,[],2);
errorbar(X,Y,Y_err,Y_err,X_err,X_err,'x')
title('Error in remaining fraction')
xlabel('Outcoupled fraction')
ylabel('Frac error in remaining pop')








