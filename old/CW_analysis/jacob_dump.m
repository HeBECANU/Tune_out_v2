%% Assuming good signals, let's set up a fit.
signal = trim_data;
fit_region = [3.6286,3.62875]*1e8;
mid_freq = median(fit_region);
fit_mask = find((signal(:,1)>fit_region(1))&(signal(:,1)<fit_region(2)));
fit_x = (signal(fit_mask,1)- mid_freq)/1e3; %Scales & centres region for better fit behavior
fit_y = signal(fit_mask,2);

%Tried a polyfit but couldn't get the root finder to work properly, this does fine
p_fun = @(p,x) p(4)+p(3)*x + p(2)*x.^2+p(1)*x.^3;
fit = fitnlm(fit_x,fit_y,p_fun,[0.5,-0.1,0.1,0]);
pc = fit.Coefficients.Estimate;
x = roots(pc);
x_real = x(imag(x)==0);


%Show me the money!
figure()
subplot(2,2,1)
scatter(signed_signal(:,1),signed_signal(:,2),'.');
title('Phase-corrected response')
xlabel('Frequency (MHz)')
ylabel('Modulation (kHz)')
subplot(2,2,2)
errorbar(trim_data(:,1),trim_data(:,2),trim_data(:,3))
title('trimmed data')
xlabel('Frequency (MHz)')
ylabel('Modulation (kHz)')
subplot(2,2,3)
plot(fit_x,fit_y,'.')
hold on
plot(fit_x,p_fun(fit.Coefficients.Estimate,fit_x));
xlabel(['\Delta=f - ',num2str(mid_freq),', GHz'])
ylabel('Modulation response')
title(['Zero crossing \Delta=',num2str(x_real),'GHz'])
subplot(2,2,4)
plot(p_fun(fit.Coefficients.Estimate,fit_x)-fit_y,'.')

%These errors are almost definitely wrong...
fit_errs = p_fun(fit.Coefficients.Estimate,fit_x)-fit_y;
SE = std(fit_errs);
MU = mean(fit_errs);

%Motivation: Find out where the curve crosses within 1sd of zero crossing,
%as this could be the ambiguity... Produces errors that seem too small to
%be believed.
x_pos_err = roots(pc-[0,0,0,SE]');
x_neg_err = roots(pc+[0,0,0,SE]');

format long
TOF = (2*mid_freq + x_real)*1e6 %TO frequency in blue
err_pos = 2*(x_real-x_pos_err(imag(x_pos_err)==0))*1e6
err_neg = 2*(x_real-x_neg_err(imag(x_neg_err)==0))*1e6
TOW = 299792458/TOF %TO wavelength in blue
W_err = -(299797458/TOF^2)*(abs(err_pos)+abs(err_neg))



%% Diagnosing spurious signals
    
%     f_test = 3.6286e8;
%     [~,t_idx] = min(abs(setpoints-f_test));
%     f_test = setpoints(t_idx);
%     test_shots = find(setpoints == f_test);
%     test_data = data.txy(test_shots);
%     figure()
%     for i=1:length(test_data)
%         ctr = ctr + 1;
%         temp_t= test_data{i}(:,1);
%         flux_t = histcounts(temp_t,hist_bin_edges);
%         [ft,c]= fft_tx(hist_bin_centers,flux_t');
%         f_c = 427;
%         [~,f_idx] = min(abs(ft(:,1)-427));
%         f_hw = 75;
%         [~,fw_idx] = min(abs(ft(:,1)-(f_c+f_hw)));
%         d_idx = fw_idx - f_idx;
%         f_window = [f_idx-d_idx:f_idx + d_idx];
%         signal = ft(f_window,2);     
%         plot(ft(f_window,1),signal)
%         hold on
%     end
%     title(['Response spectrum for f=',num2str(f_test),'Hz'])
%     xlabel('Frequency (Hz)')
%     ylabel('Response (kHz)')
    
    

    