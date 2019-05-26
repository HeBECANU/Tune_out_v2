function out=make_cal_model(anal_opts_cal,data)
%creates a model of what the trap freq of the calibration shots does over time
% TODO:
% - [ ] make sure the time is using the probe time not the labview trig time


%create a mask of shots that are calibrations and are good fits
temp_cal=data.mcp_tdc.probe.calibration';
temp_cal(isnan(temp_cal))=0;
cal_dat_mask=data.osc_fit.ok.rmse & temp_cal & ~isnan(data.wm_log.proc.probe.freq.act.mean');

temp_cal=data.mcp_tdc.probe.calibration';
temp_cal(isnan(temp_cal))=1;
probe_dat_mask=data.osc_fit.ok.rmse & ~temp_cal & ~isnan(data.wm_log.proc.probe.freq.act.mean');

%create a time vector
time_cal=data.mcp_tdc.time_create_write(cal_dat_mask,2);
time_start_cal=time_cal(1); %shift the times to start at zero for better interp perfromance
time_cal=time_cal-time_start_cal;
%create a vector of the reconstructed calibration frequencies
trap_freq_cal=data.osc_fit.trap_freq_recons(cal_dat_mask)';

%interpolate the input data to subsample by a factor of 100
time_samp_interp=linspace(min(time_cal),max(time_cal),size(time_cal,1)*1e2);
trap_freq_interp_raw=col_vec(interp1(time_cal,trap_freq_cal,time_samp_interp,'linear'));

%smooth the interpolated data
%BMH 20190523 change - moved to using gausfilt
yinterp_smooth=gaussfilt(time_samp_interp,trap_freq_interp_raw,anal_opts_cal.smooth_time);

%create a model based on interpolating the smoothed interpolated data
out.freq_drift_model=@(x) interp1(time_samp_interp,yinterp_smooth,x-time_start_cal,'linear');
out.num_shots=sum(cal_dat_mask);
out.cal_mask=cal_dat_mask;

%calculate the residuals between the model and the calibration data
model_resid=out.freq_drift_model(time_cal+time_start_cal)-trap_freq_cal;

if mean(model_resid)>std(model_resid)
    error('error mean of residuals is not within 1sd')
end
    
out.unc=std(model_resid);
mean_cal_shot_unc=mean(col_vec(data.osc_fit.trap_freq_recons_unc(cal_dat_mask)));

fprintf('%s:residual std %f vs model mean uncert %f \n',...
     mfilename,std(model_resid),mean_cal_shot_unc)


if anal_opts_cal.plot
    hour_in_s=60*60;
    x_samp=linspace(min(time_cal),max(time_cal),size(time_cal,1)*1e2);
    stfig('osc cal model','add_stack',1);
    clf
    subplot(3,1,1)
    errorbar(data.mcp_tdc.shot_num(cal_dat_mask),...
        data.osc_fit.trap_freq_recons(cal_dat_mask),data.osc_fit.trap_freq_recons_unc(cal_dat_mask)...
        ,'.k-','capsize',0,'MarkerSize',10,'linewidth',1.0)%,'r.',
    hold on
        errorbar(data.mcp_tdc.shot_num(probe_dat_mask),...
        data.osc_fit.trap_freq_recons(probe_dat_mask),data.osc_fit.trap_freq_recons_unc(probe_dat_mask)...
        ,'.b-','capsize',0,'MarkerSize',10,'linewidth',1.0)%,'r.',
    hold off
    legend('calibration data','probe data')
    xlabel('Shot Number')
    title('Input Data')
    ylabel('Measured Trap Freq (Hz)')
    
    %plot the calibration model
    subplot(3,1,2)
    plot(time_cal/hour_in_s,trap_freq_cal,'x','MarkerSize',20)
    hold on
    plot(time_samp_interp/hour_in_s,trap_freq_interp_raw,'m')
    plot(x_samp/hour_in_s,out.freq_drift_model(x_samp+time_start_cal),'k')
    legend('data','interp data','calibration model' )
    hold off
    title('Trap Freq Calibration Model')
    xlabel('experiment time (h)')
    ylabel('no probe trap freq')
    %set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1600, 900])
    %plot the calibration model
    
    subplot(3,1,3)
    errorbar(time_cal/hour_in_s,...
       model_resid,data.osc_fit.trap_freq_recons_unc(cal_dat_mask)...
        ,'.k-','capsize',0,'MarkerSize',10,'linewidth',1.0)%,'r.',
    xl=[min(time_cal),max(time_cal)]./hour_in_s;
    line(xl,[1,1]*mean(model_resid),'color','g')
    line(xl,[1,1]*mean(model_resid)+std(model_resid),'color','r')
    line(xl,[1,1]*mean(model_resid)-std(model_resid),'color','r')
    xlabel('experiment time (h)')
    ylabel('residuals')
    %set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1600, 900])
end

end


% code for ploting atom number    
%     subplot(2,1,2)
%     num_t =data.mcp_tdc.time_create_write(cal_dat_mask,2)-data.mcp_tdc.time_create_write(1,2);
%     plot(num_t/(60*60),data.mcp_tdc.num_counts(cal_dat_mask))
%     saveas(gcf,[anal_opts_cal.global.out_dir,plot_name,'.png'])
%     saveas(gcf,[anal_opts_cal.global.out_dir,plot_name,'.fig'])
%     title('Hit count trend')
%     xlabel('Time (h)')
%     ylabel('counts')
