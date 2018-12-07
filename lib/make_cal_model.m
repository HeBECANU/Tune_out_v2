function out=make_cal_model(anal_opts_cal,data)


%create a mask of shots that are calibrations and are good fits
temp_cal=data.mcp_tdc.probe.calibration';
temp_cal(isnan(temp_cal))=0;
cal_dat_mask=data.osc_fit.ok.rmse & temp_cal & ~isnan(data.wm_log.proc.probe.freq.act.mean');

%create a time vector
x_tmp=data.mcp_tdc.time_create_write(cal_dat_mask,1);
time_start_cal=x_tmp(1);
x_tmp=x_tmp-time_start_cal;
%create a vector of the reconstructed calibration frequencies
y_tmp=data.osc_fit.trap_freq_recons(cal_dat_mask)';
%interpolate the input data
xinterp=linspace(min(x_tmp),max(x_tmp),size(x_tmp,1)*1e2);
yinterp_raw=interp1(x_tmp,y_tmp,xinterp,'linear');
dx_interp=xinterp(2)-xinterp(1);
%smooth the interpolated data
kernel = gausswin(ceil(3*anal_opts_cal.smooth_time/(dx_interp)),3);
kernel=kernel/sum(kernel);
yinterp_smooth = nanconv(yinterp_raw,kernel,'edge','1d')';
%create a model based on interpolating the smoothed interpolated data
out.freq_drift_model=@(x) interp1(xinterp,yinterp_smooth,x-time_start_cal,'linear');
out.num_shots=sum(cal_dat_mask);
out.cal_mask=cal_dat_mask;
out.unc=mean(abs(out.freq_drift_model(x_tmp+time_start_cal)-y_tmp));


if anal_opts_cal.plot
    
    hour_in_s=60*60;
    x_samp=linspace(min(x_tmp),max(x_tmp),size(x_tmp,1)*1e2);
    figure(61)
    
    subplot(2,1,1)
    clf
    set(gcf,'color','w')
    plot(x_tmp/hour_in_s,y_tmp,'x','MarkerSize',20)
    hold on
    plot(xinterp/hour_in_s,yinterp_raw,'m')
    plot(x_samp/hour_in_s,out.freq_drift_model(x_samp+time_start_cal),'k')
    legend('data','interp data','calibration model' )
    hold off
    title('Trap Freq Calibration')
    xlabel('experiment time (h)')
    ylabel('no probe trap freq')
    set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1600, 900])
    plot_name='calibration_model';
    
    subplot(2,1,2)
    num_t =data.mcp_tdc.time_create_write(cal_dat_mask,2)-data.mcp_tdc.time_create_write(1,2);
    plot(num_t/(60*60),data.mcp_tdc.num_counts(cal_dat_mask))
    saveas(gcf,[anal_opts_cal.global.out_dir,plot_name,'.png'])
    saveas(gcf,[anal_opts_cal.global.out_dir,plot_name,'.fig'])
    title('Hit count trend')
    xlabel('Time (h)')
    ylabel('counts')
end

end

