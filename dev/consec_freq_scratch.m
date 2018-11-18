% figure()
% histogram(cosec_freq_dif)

%compute the difference in the same shot
in_shot_diff = cell(10);

freq_tol = 0.35;

%     freq_A = fitted_freqs{1};
%     freq_B = fitted_freqs{2};
    temp = fit_coefs{1}+1*(1/anal_opts.atom_laser.pulsedt);
    freq_A = temp(:,2,1);
    temp = fit_coefs{2}+1*(1/anal_opts.atom_laser.pulsedt);
    freq_B = temp(:,2,1);
    trap_freq_median = nanmedian(freq_A);
    mask = abs((freq_A-trap_freq_median))>freq_tol;
    mask = double(mask);
    mask(mask==1) = nan;
    freq_A = freq_A-mask;
    
    trap_freq_median = nanmedian(freq_B);
    mask = abs((freq_B-trap_freq_median))>freq_tol;
    mask = double(mask);
    mask(mask==1) = nan;
    freq_B = freq_B-mask;
    
    fprintf('\n first %u pulses \n',ii*10)
    avg_freq_A = nanmean(freq_A);
    avg_freq_B = nanmean(freq_B);
    in_shot_diff{ii} = (freq_A).^2-(freq_B-avg_freq_B+avg_freq_A).^2;
    fprintf('comparison std: %f \n',nanstd((freq_A).^2-(freq_B).^2))
    figure
    histogram((freq_A).^2-(freq_B).^2,100)
    xlabel('In shot square difference')
    ylabel('count')
    %number of faled fitts
    fprintf('failed shots: %u \n', sum(isnan(in_shot_diff{ii})))
    display(strcat([num2str(nanmean(in_shot_diff{ii})),' +/- ',num2str(nanstd(in_shot_diff{ii}))]))


%display(strcat('Consecutive shot difference= ',num2str(mean(abs(cosec_freq_dif))), char(177) ,num2str(std(cosec_freq_dif))))
%%

    display(strcat(['std A = ',num2str(nanstd(fitted_freqs{1}.^2)),', std B = ',num2str(nanstd(fitted_freqs{2}.^2))]))


%%
temp = ones(147,3);
for ii = 1:3
    for jj = 1:147
        temp(jj,ii) = fitted_freqs{1,ii}(jj);
    end
end

figure
plot(temp')

%%
figure
drift_A=freq_drift{1};
drift_B=freq_drift{2};
plot(drift_A)
hold on
plot(drift_B)

%% compare fitted frequency over time for different sample sizes
figure
r = 19;
fit_param_vals = fit_coefs{2};
plot(t_var,fit_param_vals(r,2,1).*exp(-t_var.*fit_param_vals(r,7,1)))
hold on
fit_param_vals = fit_coefs{3};
plot(t_var,fit_param_vals(r,2,1).*exp(-t_var.*fit_param_vals(r,7,1)))
fit_param_vals = fit_coefs{1};
plot(t_var,fit_param_vals(r,2,1).*exp(-t_var.*fit_param_vals(r,7,1)))

%%
