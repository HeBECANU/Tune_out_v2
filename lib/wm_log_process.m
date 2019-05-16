function wm_log_processed=wm_log_process(anal_opts,data)
%this is so quick it barely needs cacheing
cache_opts=[];
cache_opts.verbose=0;
%cache_opts.force_cache_load=wm_log_import_opts.force_load_save;
%wm_log_import_opts=rmfield(wm_log_import_opts,'force_load_save');
cache_opts.force_recalc=anal_opts.wm_log.force_reimport;
cache_opts.save_compressed=true;
anal_opts.wm_log=rmfield(anal_opts.wm_log,'force_reimport');
%create a temp copy of data to limit the effect that changes have on the cache function
%this is the limitation of the pass data structure approach, the scope can quickly become to large
sub_data=[];
sub_data.wm_log.raw=data.wm_log.raw;
sub_data.mcp_tdc.time_create_write=data.mcp_tdc.time_create_write;
sub_data.labview.time=data.labview.time;

outputs=function_cache(cache_opts,@wm_log_process_core,{anal_opts,sub_data});
wm_log_processed=outputs{1};
end

function wm_log_processed=wm_log_process_core(anal_opts,data)
%mean_shot_duration=mean(diff(data.labview.time(:)));

wm_log_processed=[];
fprintf('finding mean wavelengths for files %04u:%04u',size(data.mcp_tdc.time_create_write(:,1),1),0)
data.wm_log.raw.feedback.posix_time = sort(data.wm_log.raw.feedback.posix_time);
data.wm_log.raw.read_all_adc.posix_time = sort(data.wm_log.raw.read_all_adc.posix_time);
if ~issorted(data.wm_log.raw.feedback.posix_time) || ~issorted(data.wm_log.raw.read_all_adc.posix_time)
    error('binary search wont work, keep it ordered plz')
end
iimax=size(data.mcp_tdc.time_create_write(:,1),1);

wm_log_processed.probe.freq.set=nan(iimax,1);
wm_log_processed.probe.freq.act.mean=nan(iimax,1);
wm_log_processed.probe.freq.act.std=nan(iimax,1);
wm_log_processed.probe.freq.rvb.mean=nan(iimax,1);

%GUILTY UNTILL PROVEN INOCENT  !!!!!!!!!!!!!!!!!
%data.mcp_tdc.probe.ok.freq=false(iimax,1); %frequency reading
%data.mcp_tdc.probe.ok.rvb=false(iimax,1);  %2r-b check
%data.mcp_tdc.probe.ok.ecd_pd=false(iimax,1);  %ecd pd value

wm_log_processed.ok.freq=false(iimax,1); %frequency reading
wm_log_processed.ok.rvb=false(iimax,1);  %2r-b check
wm_log_processed.ok.ecd_pd=false(iimax,1);  %ecd pd value

sfigure(11);
set(gcf,'color','w')

imax=min([size(data.labview.time,2),size(data.mcp_tdc.time_create_write,1)]);
time_diff=data.mcp_tdc.time_create_write(1:imax,2)'-anal_opts.dld_aquire-anal_opts.trig_dld-...
        data.labview.time(1:imax);
mean_delay_labview_tdc=mean(time_diff);

for ii=1:iimax
    est_labview_start=data.mcp_tdc.time_create_write(ii,2)...
            -anal_opts.trig_dld-anal_opts.dld_aquire-mean_delay_labview_tdc;
        
    %find the labview update that was nearest to the tdc_time minus the tdc_trig time
    %this should be the labview itteration that made this data 
    %is this the best way to do it?
    [time_nearest_lv,idx_nearest_lv]=closest_value(data.labview.time,est_labview_start);
    time_diff=data.labview.time(idx_nearest_lv)-est_labview_start;
    
    abs_anal_opts.trig_dld_on_main_comp=data.labview.time(idx_nearest_lv)+anal_opts.trig_dld;
    time_lower=abs_anal_opts.trig_dld_on_main_comp+anal_opts.atom_laser.t0-anal_opts.global.fall_time-0.5; %the probe turns on
    time_upper=abs_anal_opts.trig_dld_on_main_comp+anal_opts.atom_laser.t0-anal_opts.global.fall_time+anal_opts.wm_log.time_probe+0.5; %when the probe beam goes off
    
    %CHECK PD VALUE
    %check if there were some readings of the ecd output voltage withing time_padding seconds of this period
    %the padding is because EDC readings are only taken every 3s
    ecd_pd_mask_idx=fast_sorted_mask(data.wm_log.raw.read_all_adc.posix_time,...
                time_lower-anal_opts.wm_log.time_pd_padding,...
                time_upper+anal_opts.wm_log.time_pd_padding);
    ecd_pd_measurments=max(0,ecd_pd_mask_idx(2)-ecd_pd_mask_idx(1));
    pd_readings=data.wm_log.raw.read_all_adc.value10(ecd_pd_mask_idx(1):ecd_pd_mask_idx(2));  
    if sum(ecd_pd_measurments)>=2 && sum(pd_readings<anal_opts.wm_log.ecd_volt_thresh)==0 %if there are not a few during this time then throw an error
        wm_log_processed.ok.ecd_pd(ii)=true;
    end
    
    %CHECK wm freq
    %could check if close to set pt, but for now will just check if flat
    wm_red_mask_idx=fast_sorted_mask(data.wm_log.raw.feedback.posix_time,time_lower,time_upper);
    red_freqs=data.wm_log.raw.feedback.actual(wm_red_mask_idx(1):wm_red_mask_idx(2));
    wm_log_processed.probe.freq.set(ii)=mean(data.wm_log.raw.feedback.setpt(wm_red_mask_idx(1):wm_red_mask_idx(2)));
    wm_log_processed.probe.freq.act.mean(ii)=mean(red_freqs);
    wm_log_processed.probe.freq.act.std(ii)=std(red_freqs);
    wm_log_processed.probe.freq.error(ii)=false;
    %BMH todo: if there 
    if numel(red_freqs)>5 && ... %need more than 5 measurments to say anything about stability
            std(red_freqs)<anal_opts.wm_log.red_sd_thresh && range(red_freqs)<anal_opts.wm_log.red_range_thresh
       wm_log_processed.ok.freq(ii)=true;
    end
    
    %check that the blue is ~ 2x the red freq
    wm_blue_mask_idx=fast_sorted_mask(data.wm_log.raw.blue_freq.posix_time,...
                time_lower-anal_opts.wm_log.time_blue_padding,...
                time_upper+anal_opts.wm_log.time_blue_padding);
    blue_freqs=data.wm_log.raw.blue_freq.value(wm_blue_mask_idx(1):wm_blue_mask_idx(2));
    
    wm_log_processed.probe.freq.rvb.mean(ii)=mean(blue_freqs)-wm_log_processed.probe.freq.act.mean(ii)*2;
    if abs(wm_log_processed.probe.freq.rvb.mean(ii))<anal_opts.wm_log.rvb_thresh
        wm_log_processed.ok.rvb(ii)=true;
    end
    %do some plots if plot_all or if plot_failed what failed
    if anal_opts.wm_log.plot_all || (anal_opts.wm_log.plot_failed &&...
       (~wm_log_processed.ok.freq(ii) || ~wm_log_processed.ok.rvb(ii) ||...
       ~wm_log_processed.ok.ecd_pd(ii)))
        
        sfigure(11);
        subplot(3,1,1)
        plot(data.wm_log.raw.read_all_adc.posix_time(ecd_pd_mask_idx(1):...
            ecd_pd_mask_idx(2))+anal_opts.atom_laser.t0-...
            abs_anal_opts.trig_dld_on_main_comp,pd_readings)
        xlabel('Probe on time(s)')
        ylabel('Doubler PD voltage')
        title('Doubler output check')
        subplot(3,1,2)
        tmp_wav_avg=mean(red_freqs);
        tmp_wav_cen=red_freqs-tmp_wav_avg;
        plot(anal_opts.atom_laser.t0+ data.wm_log.raw.feedback.posix_time(wm_red_mask_idx(1):wm_red_mask_idx(2))...
            -abs_anal_opts.trig_dld_on_main_comp,tmp_wav_cen)
        title(sprintf('Red freq - %.1fMHz',tmp_wav_avg))
        xlabel('Probe on time(s)')
        ylabel('Red Freq (MHz)')
        yl=ylim;

        
        subplot(3,1,3)
        %should really interpolate here to give the difference properly as
        %a function of time
        title('2red(mean)-blue difference')
        tmp_rvb=blue_freqs-wm_log_processed.probe.freq.act.mean(ii)*2;
        plot(anal_opts.atom_laser.t0+ data.wm_log.raw.blue_freq.posix_time(wm_blue_mask_idx(1):wm_blue_mask_idx(2))...
            -abs_anal_opts.trig_dld_on_main_comp,tmp_rvb)
        xlabel('Probe on time(s)')
        ylabel('2r-b Freq (MHz)')
        title('Blue red difference')
        pause(1e-6)
    end
    if mod(ii,1e1)==0,fprintf('\b\b\b\b%04u',ii),end
end
fprintf('...Done\n')




end