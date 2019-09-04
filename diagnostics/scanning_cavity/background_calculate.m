function out_struct=background_calculate(in_struct)
% a lot of code taken from is_laser_single_mode 
% https://github.com/brycehenson/is_laser_single_mode
% can only handle sawtooth with +ve polarity currently

cli_format_text('background_calculate','h',2)

%------------- BEGIN USER VAR --------------
in_struct.pzt_scan_period=1;
sawtooth_polarity_positive=true;
st_pt_min_pzt_factor=0.5;
estimated_finesse=100;
%------------- END USER VAR --------------

%------------- BEGIN CODE --------------

out_struct=[];
% %save('_T__before_background_calc.mat')

%% measure the scan freq
% not really needed at this point

if isfield(in_struct,'pzt_scan_period')
    tmax=min(in_struct.times(1)+in_struct.pzt_scan_period*10,in_struct.times(end));
    [~,idx_max]=closest_value(in_struct.times,tmax); %or could use fast_sorted_mask
    pztfreq_subset_data=[in_struct.times(1:idx_max),in_struct.pzt_voltage(1:idx_max)];
    %clear('tmax','idx_max')
else %otherwise we do the fft on the whole thing
    pztfreq_subset_data=[in_struct.times,in_struct.pzt_voltage];
end
% 
% 
% 
% df_args=[];
% df_args.components_min_amp=1e-2;
% df_args.freq_limits=[0.1,10].*(1/in_struct.pzt_scan_period);
% dom_freq_out=dominant_freq_components(pztfreq_subset_data(:,1),pztfreq_subset_data(:,2),[])

%%
pzt_scan_std=std(pztfreq_subset_data(:,2));
pzt_scan_range=range(pztfreq_subset_data(:,2));

pzt_scan_period=in_struct.pzt_scan_period;

time_range=in_struct.times(end)- in_struct.times(1);
scans_in_data=time_range/pzt_scan_period;

%% get the stationary points of the pzt waveform

cli_format_text('splitting up by pzt scans','h',3)

st_pts=stpt_and_level_xing(in_struct.times,in_struct.pzt_voltage,...
                                        in_struct.pzt_filt_factor_deriv*pzt_scan_period,nan);

%st_pts.value=in_struct.pzt_voltage(st_pts.idx);
pzt_diff_between_st_pts=diff(st_pts.xval);
%select pairs that have a large enough difference
scan_mask=pzt_diff_between_st_pts>pzt_scan_range*st_pt_min_pzt_factor;
if strcmp(in_struct.scan_type,'sawtooth') 
    %select the pairs of stationary points that give a differenc that is the right sign(positive if sawtooth_polarity_positive==true)
    if sawtooth_polarity_positive
        scan_mask=scan_mask & pzt_diff_between_st_pts>0;
    else
        scan_mask=scan_mask & pzt_diff_between_st_pts<0;
    end
end
scan_stpt_idxs=find(scan_mask); %find the first scan that satisfies the conditions
%make a matrix with start and stop indicies for all scans (for the long slope)
pzt_scans=[];
pzt_scans.start.idx=st_pts.idx(scan_stpt_idxs);
pzt_scans.start.time=in_struct.times(pzt_scans.start.idx);
pzt_scans.stop.idx=st_pts.idx(scan_stpt_idxs+1);                                   
pzt_scans.stop.time=in_struct.times(pzt_scans.stop.idx);


%% find all the pd peaks
pd_thresh_fac=0.1;

cli_format_text('finding all pd peaks (rough)','h',3)
pd_dat_uncmp=in_struct.pd_voltage(:,1);
pd_dat_cmp=in_struct.pd_voltage(:,2);
pd_dat_uncmp_filt=gaussfilt(in_struct.times,pd_dat_uncmp,...
    in_struct.pd_filt_factor*in_struct.pzt_scan_period);

pd_uncmp_range=range(pd_dat_uncmp_filt);



pzt_voltage=in_struct.pzt_voltage;
ptz_voltage_filt=gaussfilt(in_struct.times,pzt_voltage,...
    in_struct.ptz_filt_factor_pks*in_struct.pzt_scan_period);


[pks_pd,pd_pks.idx] = findpeaks(pd_dat_uncmp_filt,'MinPeakHeight', pd_thresh_fac*pd_uncmp_range);
pd_pks.pd=pks_pd;
pd_pks.time=in_struct.times(pd_pks.idx);
pd_pks.pzt=ptz_voltage_filt(pd_pks.idx);

pzt_scans.start.pzt=ptz_voltage_filt(pzt_scans.start.idx);
pzt_scans.stop.pzt=ptz_voltage_filt(pzt_scans.stop.idx);


%%
stfig('peak finding diagnostic','add_stack',1);
clf
plot(in_struct.times,ptz_voltage_filt,'b')
hold on
plot(in_struct.times,pd_dat_uncmp_filt,'r')
plot(pd_pks.time,pd_pks.pd,'xr')
plot(pzt_scans.start.time,pzt_scans.start.pzt,'xg')
plot(pzt_scans.stop.time,pzt_scans.stop.pzt,'xk')
hold off


pzt_scans.start.time=in_struct.times(pzt_scans.start.idx);
pzt_scans.stop.idx=st_pts.idx(scan_stpt_idxs+1);                                   
pzt_scans.stop.time=in_struct.times(pzt_scans.stop.idx);

%%
cli_format_text('...and he separated the light from the darkness','h',3) %

%lets try to figure out the bounds on the regions where we should calculate the light on & light off information

pd_peak_time_diff=diff(pd_pks.time);
pd_peak_time_diff_trimmed=pd_peak_time_diff(pd_peak_time_diff<(pzt_scan_period/2));
stfig('pd peak time differences','add_stack',1);
% calculate the trimed mean of the time differences
pd_peak_time_diff_mean=mean(pd_peak_time_diff_trimmed);
pd_peak_time_diff_std=std(pd_peak_time_diff_trimmed);



%now use this to split up the light & dark

%find when there is a large gap in times between peaks
%note:this is not perfect bc if the blanking of the pd amplifier is too long it wont work
% however here its fine
pd_peak_gap_logic=pd_peak_time_diff>pd_peak_time_diff_mean+pd_peak_time_diff_std*3;

is_light_on=consecutive_true(~pd_peak_gap_logic,3);
light_on_sections=[];
light_on_sections.pk_idx.on=strfind(is_light_on,[0,1])+1;
if is_light_on(1), light_on_sections.pk_idx.on=[1,light_on_sections.pk_idx.on]; end

light_on_sections.pk_idx.off=strfind(is_light_on,[1,0]);
if is_light_on(end), light_on_sections.pk_idx.off=[light_on_sections.pk_idx.off,numel(is_light_on)]; end

% we can now turn this into on and off transtion times
light_on_sections.time.start=pd_pks.time(light_on_sections.pk_idx.on);
light_on_sections.time.stop=pd_pks.time(light_on_sections.pk_idx.off);
% for the light on & light off interogation times we need to padd things differently

%trim off the first and last peak just to be extra safe

light_on_sections.time_padd.start=light_on_sections.time.start+20*pd_peak_time_diff_mean/estimated_finesse;
light_on_sections.time_padd.stop=light_on_sections.time.stop-20*pd_peak_time_diff_mean/estimated_finesse;

% convert to mask
light_on_sections.mask=false*in_struct.times;
for scan_idx=1:numel(light_on_sections.time_padd.start)
    idxs=fast_sorted_mask(in_struct.times,...
        light_on_sections.time_padd.start(scan_idx),...
        light_on_sections.time_padd.stop(scan_idx));
    light_on_sections.mask(idxs(1):idxs(2))=true;
end




%% find the dark sections

pd_peak_gap_logic=pd_peak_time_diff>pd_peak_time_diff_mean*2;
light_off_sections=[];
light_off_sections.pk_idx.start=find(pd_peak_gap_logic);
light_off_sections.pk_idx.stop=light_off_sections.pk_idx.start+1;

light_off_sections.time.start=pd_pks.time(light_off_sections.pk_idx.start);
light_off_sections.time.stop=pd_pks.time(light_off_sections.pk_idx.stop);

% deal with if the file ends at a light off condition
if pd_pks.time(end)+pd_peak_time_diff_mean*4<in_struct.times(end)
    light_off_sections.time.start=[light_off_sections.time.start,pd_pks.time(end)];
    light_off_sections.time.stop=[light_off_sections.time.stop,in_struct.times(end)];
end

light_off_sections.time_padd.start=light_off_sections.time.start+2*pd_peak_time_diff_mean;
light_off_sections.time_padd.stop=light_off_sections.time.stop-2*pd_peak_time_diff_mean;

% converting to masks
% to make calculations easier i will convert these sections into a mask
light_off_sections.mask=false*in_struct.times;
for scan_idx=1:numel(light_off_sections.time_padd.start)
    idxs=fast_sorted_mask(in_struct.times,...
        light_off_sections.time_padd.start(scan_idx),...
        light_off_sections.time_padd.stop(scan_idx));
    light_off_sections.mask(idxs(1):idxs(2))=true;
end

%% now find the mean(and std) of the light off voltage

% first make a mask of the scans
pzt_scans.mask=false*in_struct.times;
for scan_idx=1:numel(pzt_scans.start.idx)
    pzt_scans.mask(pzt_scans.start.idx(scan_idx):pzt_scans.stop.idx(scan_idx))=true;
end
%% Diagnostic plot
stfig('background selection diagnostic','add_stack',1);
clf

pd_scaling_fun=@(x) 3*log10(1+log10(x+1));
plot(in_struct.times,pd_scaling_fun(pd_dat_uncmp-min(pd_dat_uncmp)))
hold on
plot(in_struct.times,pzt_voltage*0.07)
plot(in_struct.times,light_off_sections.mask)
plot(in_struct.times,pzt_scans.mask)
hold off




%%


mask_light_off_during_scan=light_off_sections.mask &  pzt_scans.mask;
pd_vals_uncmp_light_off=pd_dat_uncmp(mask_light_off_during_scan);
times_light_off=in_struct.times(mask_light_off_during_scan);

%%
stfig('background drift diagnostic','add_stack',1);
clf
%plot(times_light_off,pd_vals_uncmp_light_off)

% gaussfilt is too slow for non ordered data
%pd_vals_uncmp_light_off_smooth=gaussfilt(times_light_off,pd_vals_uncmp_light_off,1);
pd_vals_uncmp_light_off_smooth=movmean(pd_vals_uncmp_light_off,1e4);
plot(times_light_off,pd_vals_uncmp_light_off_smooth)

%in future i would like to make a smoothing model so that the background is not assumed to be uniform



%%

out_struct.pd_light_off=[];
out_struct.pd_light_off.num_samp=numel(pd_vals_uncmp_light_off);

out_struct.pd_light_off.uncmp.mean=mean(pd_vals_uncmp_light_off);
out_struct.pd_light_off.uncmp.std=std(pd_vals_uncmp_light_off);
out_struct.pd_light_off.uncmp.ste=out_struct.pd_light_off.uncmp.std/sqrt(out_struct.pd_light_off.num_samp);

pd_vals_cmp_light_off=pd_dat_cmp(mask_light_off_during_scan);
out_struct.pd_light_off.cmp.mean=mean(pd_vals_cmp_light_off);
out_struct.pd_light_off.cmp.std=std(pd_vals_cmp_light_off);
out_struct.pd_light_off.cmp.ste=out_struct.pd_light_off.cmp.std/sqrt(out_struct.pd_light_off.num_samp);

%% compare the histograms to see if it is worth using the compressed channel
do_no_light_cmp_comparison=0;
if do_no_light_cmp_comparison

    smoothing=1e-6
    sh_out=smooth_hist(pd_vals_uncmp_light_off,'sigma',smoothing);

    stfig('voltage density')
    clf
    plot( sh_out.bin.centers,sh_out.count_rate.smooth/max(sh_out.count_rate.smooth))
    hold on
    sh_out=smooth_hist(pd_vals_cmp_light_off,'sigma',smoothing);
    plot( sh_out.bin.centers,sh_out.count_rate.smooth/max(sh_out.count_rate.smooth))

end
% ok so they seem very similar

% now we are done calculating the background

%% light on measurment
% take all the selected peaks

peaks_mask=light_on_sections.mask(pd_pks.idx);
pd_peak_pairs.idx.start=pd_pks.idx(logical(peaks_mask));
pd_peak_pairs.idx.start=pd_peak_pairs.idx.start(1:end-1);
pd_peak_pairs.idx.stop=pd_pks.idx(logical(peaks_mask));
pd_peak_pairs.idx.stop=pd_peak_pairs.idx.stop(2:end);

pd_peak_pairs.time.start=in_struct.times(pd_peak_pairs.idx.start);
pd_peak_pairs.time.stop=in_struct.times(pd_peak_pairs.idx.stop);


% now we go through and veto the pair if
% - there is a scan end in between
% - the light goes off in between
pd_peak_pairs_mask=false(numel(pd_peak_pairs.idx.start),1);
for ii=1:numel(pd_peak_pairs.idx.start)
    is_light_on_between=sum(pd_peak_pairs.time.start(ii)<light_on_sections.time_padd.stop &... 
                    pd_peak_pairs.time.stop(ii)>light_on_sections.time_padd.stop)==0;
    pzt_scan_continous=sum(pd_peak_pairs.time.start(ii)<pzt_scans.stop.time &... 
                    pd_peak_pairs.time.stop(ii)>pzt_scans.stop.time)==0;
    pzt_scan_continous=pzt_scan_continous & sum(pd_peak_pairs.time.start(ii)<pzt_scans.start.time &... 
                    pd_peak_pairs.time.stop(ii)>pzt_scans.stop.time)==0;
    pd_peak_pairs_mask(ii)=is_light_on_between & pzt_scan_continous;
    
end

pd_peak_pairs.idx.start=pd_peak_pairs.idx.start(pd_peak_pairs_mask);
pd_peak_pairs.idx.stop=pd_peak_pairs.idx.stop(pd_peak_pairs_mask);
pd_peak_pairs.time.start=pd_peak_pairs.time.start(pd_peak_pairs_mask);
pd_peak_pairs.time.stop=pd_peak_pairs.time.stop(pd_peak_pairs_mask);
pd_peak_pairs.pd_dat_uncmp.start=pd_dat_uncmp_filt(pd_peak_pairs.idx.start);
pd_peak_pairs.pd_dat_uncmp.stop=pd_dat_uncmp_filt(pd_peak_pairs.idx.stop);

%% Diagnostic plot
stfig('peak selection diagnostic','add_stack',1);
clf
min_pd_val=min(pd_dat_uncmp_filt)
plot(in_struct.times,log10(pd_dat_uncmp_filt-min_pd_val+1))
hold on
plot(in_struct.times,pzt_voltage*0.1)
plot(in_struct.times,light_off_sections.mask)
plot(in_struct.times,pzt_scans.mask)
peak_times_mat=cat(2,pd_peak_pairs.time.start,pd_peak_pairs.time.stop)';
peak_pd_mat=cat(2,pd_peak_pairs.pd_dat_uncmp.start,pd_peak_pairs.pd_dat_uncmp.stop)';
plot(peak_times_mat,log10(peak_pd_mat-min_pd_val+1))
hold off


%%
int_range_cen=0.2;
int_range_peak=0.05;
int_cen_vec=nan(numel(pd_peak_pairs.idx.start),1);
int_peak_vec=int_cen_vec;
for ii=1:numel(pd_peak_pairs.idx.start)
    mean_time=mean([pd_peak_pairs.time.start(ii),pd_peak_pairs.time.stop(ii)]);
    time_diff=range([pd_peak_pairs.time.start(ii),pd_peak_pairs.time.stop(ii)]);
    
    int_time_lims=mean_time+time_diff*int_range_cen*[-1,1]*0.5;
    int_mask=fast_sorted_mask(in_struct.times,int_time_lims(1),int_time_lims(2));
    int_cen_xval=in_struct.times(int_mask(1):int_mask(2));
    int_cen_yval= pd_dat_uncmp(int_mask(1):int_mask(2))- out_struct.pd_light_off.uncmp.mean;
    int_cen_vec(ii)=trapz(int_cen_xval,int_cen_yval);
    
    
    int_time_lims=pd_peak_pairs.time.start(ii)+time_diff*int_range_peak*[-1,1]*0.5;
    int_mask=fast_sorted_mask(in_struct.times,int_time_lims(1),int_time_lims(2));
    int_peak_start_xval=in_struct.times(int_mask(1):int_mask(2));
    int_peak_start_yval= pd_dat_uncmp(int_mask(1):int_mask(2))- out_struct.pd_light_off.uncmp.mean;
    int_peak_start=trapz(int_peak_start_xval,int_peak_start_yval);
    
    
    int_time_lims=pd_peak_pairs.time.stop(ii)+time_diff*int_range_peak*[-1,1]*0.5;
    int_mask=fast_sorted_mask(in_struct.times,int_time_lims(1),int_time_lims(2));
    int_peak_stop_xval=in_struct.times(int_mask(1):int_mask(2));
    int_peak_stop_yval= pd_dat_uncmp(int_mask(1):int_mask(2))- out_struct.pd_light_off.uncmp.mean;
    if numel(int_peak_stop_xval)~=numel(int_peak_stop_yval)
        fprintf('debug')
    end
    int_peak_stop=trapz(int_peak_stop_xval,int_peak_stop_yval);
    
    int_peak_vec(ii)=mean([int_peak_start,int_peak_stop]);
    
    if in_struct.plot.all
        stfig('integeration diagnostic','add_stack',1);
        area(int_peak_stop_xval-mean_time,int_peak_stop_yval)
        hold on
        area(int_peak_start_xval-mean_time,int_peak_start_yval)
        area(int_cen_xval-mean_time,int_cen_yval)
        hold off
        title(sprintf('peak pair at %f s',mean_time))
        drawnow
        pause(0.001)
    end
    
    
    
end


%%
power_ratio=int_cen_vec./int_peak_vec



cli_format_text('done','h',3) %and he separated the light from the darkness





end