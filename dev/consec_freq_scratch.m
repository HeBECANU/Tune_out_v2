% figure()
% histogram(cosec_freq_dif)

%compute the difference in the same shot
in_shot_diff = cell(10);

freq_tol = 0.35;

for ii = 7
    freq_A = fitted_freqs{1}{ii};
    freq_B = fitted_freqs{2}{ii};
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
    histogram((freq_A).^2-(freq_B).^2)
    %number of faled fitts
    fprintf('failed shots: %u \n', sum(isnan(in_shot_diff{ii})))
    display(strcat([num2str(nanmean(in_shot_diff{ii})),' +/- ',num2str(nanstd(in_shot_diff{ii}))]))
end

%display(strcat('Consecutive shot difference= ',num2str(mean(abs(cosec_freq_dif))), char(177) ,num2str(std(cosec_freq_dif))))
%%
for ii=1:7
    display(strcat(['std A = ',num2str(nanstd(fitted_freqs{1}{ii}.^2)),', std B = ',num2str(nanstd(fitted_freqs{2}{ii}.^2))]))
end
