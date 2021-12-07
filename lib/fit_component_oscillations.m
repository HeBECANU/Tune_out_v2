function fit_component_oscillations(segmentd_al,options)


num_shots=size(segmentd_al.segmented_vel_xyz.high_flux.mean,1);
num_pulses=size(segmentd_al.segmented_vel_xyz.high_flux.mean,2);

pulse_num=1:num_pulses;

for shot=1:num_shots
    %%
    dim=3;
    high_mean_vel=segmentd_al.segmented_vel_xyz.high_flux.mean(shot,:,dim);
    high_ste_vel=segmentd_al.segmented_vel_xyz.high_flux.ste(shot,:,dim);
    low_mean_vel=segmentd_al.segmented_vel_xyz.low_flux.mean(shot,:,dim);
    low_ste_vel=segmentd_al.segmented_vel_xyz.low_flux.ste(shot,:,dim);
    close_low_mean_vel=segmentd_al.segmented_vel_xyz.close_low_flux.mean(shot,:,dim);
    close_low_ste_vel=segmentd_al.segmented_vel_xyz.close_low_flux.ste(shot,:,dim);
    
    %
    stfig('component osc compare');
    clf
    subplot(2,1,1)
    errorbar(segmentd_al.time_cen,low_mean_vel,low_ste_vel,'or-','Capsize',0,'Linewidth',1.5) % ,'LineStyle','none'
    hold on
    errorbar(segmentd_al.time_cen',high_mean_vel,high_ste_vel,'ok-','Capsize',0,'Linewidth',1.5) % 'LineStyle','none'
    errorbar(segmentd_al.time_cen',close_low_mean_vel,close_low_ste_vel,'og-','Capsize',0,'Linewidth',1.5) % 'LineStyle','none'
    hold off
    subplot(2,1,2)
    low_fft=fft_tx(segmentd_al.time_cen,low_mean_vel-mean(low_mean_vel));
    plot(low_fft(1,:),abs(low_fft(2,:)),'k')
    hold on
    high_fft=fft_tx(segmentd_al.time_cen,high_mean_vel-mean(high_mean_vel));
    plot(high_fft(1,:),abs(high_fft(2,:)),'r')
    close_low_fft=fft_tx(segmentd_al.time_cen,close_low_mean_vel-mean(close_low_mean_vel));
    plot(close_low_fft(1,:),abs(close_low_fft(2,:)),'g')
    hold off
    
%     [xcorr_amp,xcorr_lags]=xcorr(high_mean_vel-mean(high_mean_vel),low_mean_pos-mean(low_mean_pos));
%      plot(xcorr_lags,xcorr_amp)


end



end