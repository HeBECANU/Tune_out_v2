%checking that matlab can happily injest the data
import_opts.dir='Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\';
ai_log_import_opts.dir=import_opts.dir;
ai_log_import_opts.force_reimport=false;
wm_log_name='log_analog_in_';
wm_logs=dir([ai_log_import_opts.dir,wm_log_name,'*.txt']);
ai_log_import_opts.names={wm_logs.name};
path=strcat(ai_log_import_opts.dir,ai_log_import_opts.names{end-0});
fid = fopen(path,'r');
ai_log_file_cells=textscan(fid,'%s','Delimiter','\n');
fclose(fid);

ai_dat=jsondecode(ai_log_file_cells{1}{1});
samples= size(ai_dat.Data,2);
sr=ai_dat.sample_rate;


pd_set=2;
mean_diff_thesh=30e-3;
std_thesh=20e-3;
time_start=0.2;
time_stop=2;

sampl_start=max(1,ceil(time_start*sr));
sampl_stop=min(samples,ceil(time_stop*sr));
probe_pd_during_meas=ai_dat.Data(1,sampl_start:sampl_stop);

if abs(mean(probe_pd_during_meas)-pd_set)>mean_diff_thesh
    fprintf('probe beam pd value wrong!!!!!!!!!\n')
end
if std(probe_pd_during_meas)>std_thesh
    fprintf('probe beam pd noisy!!!!!!!!!\n')
end
sfigure(2);
clf
plot(probe_pd_during_meas)





%%

%%now thats the dat is injested lets try and see if its multimode
scan_time=25e-3;
test_times=linspace(0,3-scan_time,50);
thresh_cmp_peak=20e-3;%volts
pzt_division=0.2498;
peak_distance_thresh_cmp_full_min=0.3;%volts difference full-cmp peak in order to be considered new peak
peak_dist_min_pass=4.5;%minimum (min difference)between peaks for the laser to be considered single mode
plots=false;

test=true; 
ii=0;
iimax=numel(test_times);
sweep=[];
single_mode_vec=NaN(1,iimax);
fprintf('%03i',0)
while test
ii=ii+1;
time_start=test_times(ii);
fprintf('\b\b\b%03i',ii)

sampl_start=max(1,ceil(time_start*sr));
sampl_stop=min(samples,sampl_start+ceil(scan_time*sr));
if sampl_stop-sampl_start>=scan_time*sr %check that we have enough points to work with
sub_dat_raw=ai_dat.Data(3,sampl_start:sampl_stop)/pzt_division;
window=gausswin(30,4);
window=window/sum(window(:));%normalize
sub_dat_smooth=conv(sub_dat_raw,window,'same');
sub_dat_grad_smooth=diff(sub_dat_smooth)*sr;
pos_slope_mask=sub_dat_grad_smooth>0;


%now we want to select the first positive going sgement that is long enough

x=pos_slope_mask';
% Find the transitions
down= strfind(x',[1 0]);
ups= strfind(x',[0 1]);

%we take the ramp to investigate as the between the first rising edge of the mask and the next
%falling edge

sweep{ii}.start_idx=ups(1);
%find the next falling edge
sweep{ii}.stop_idx=down(find(down>sweep{ii}.start_idx,1)); %fuck so elegant, fuck this is so much beter than LabView

sweep{ii}.pzt_raw=sub_dat_raw(sweep{ii}.start_idx:sweep{ii}.stop_idx);
sweep{ii}.pzt_smooth=sub_dat_smooth(sweep{ii}.start_idx:sweep{ii}.stop_idx);
sweep{ii}.pd_full_raw=ai_dat.Data(2,sampl_start+sweep{ii}.start_idx:sampl_start+sweep{ii}.stop_idx);
sweep{ii}.pd_cmp_raw=ai_dat.Data(4,sampl_start+sweep{ii}.start_idx:sampl_start+sweep{ii}.stop_idx);
sweep{ii}.pd_cmp_raw=sweep{ii}.pd_cmp_raw-median(sweep{ii}.pd_cmp_raw);
%threshold out all data that is 1/10 the max value
thesh=0.07*max(sweep{ii}.pd_full_raw);
pd_mask=sweep{ii}.pd_full_raw>thesh;

cmp_multiplier_disp=50; %multiplier to display the compressed data better


[pks_val,locs] = findpeaks(sweep{ii}.pd_full_raw,'MinPeakHeight',thesh);
sweep{ii}.pks.full.pd=pks_val;
sweep{ii}.pks.full.pzt=sweep{ii}.pzt_smooth(locs);

[pks_val,locs] = findpeaks(sweep{ii}.pd_cmp_raw,'MinPeakHeight',thresh_cmp_peak);
sweep{ii}.pks.cmp.pd=pks_val;
sweep{ii}.pks.cmp.pzt=sweep{ii}.pzt_smooth(locs);


%find the distance in pzt voltage to the nearest full peak
min_dist_cmp_to_full=arrayfun(@(x) min(abs(x-sweep{ii}.pks.full.pzt)),sweep{ii}.pks.cmp.pzt);
%peak_distance_thresh_cmp_full_min
sweep{ii}.pks.cmp.pzt=sweep{ii}.pks.cmp.pzt(min_dist_cmp_to_full>peak_distance_thresh_cmp_full_min);
sweep{ii}.pks.cmp.pd=sweep{ii}.pks.cmp.pd(min_dist_cmp_to_full>peak_distance_thresh_cmp_full_min);




if plots
    sfigure(1);
    clf
    subplot(2,2,1)
    title('smoothing pzt v')
    set(gcf,'color','w')
    time=(1:numel(sub_dat_raw))/sr;
    plot(time*1e3,sub_dat_raw,'r')
    hold on
    plot(time*1e3,sub_dat_smooth,'k')
    hold off
    xlabel('time (ms)')
    ylabel('volts (v)')

    subplot(2,2,2)
    %select a single scan
    %to do this we will mask out the positive slope
    diff_sub_raw=diff(sub_dat_raw)*sr;
    diff_time=time(1:end-1)+0.5*(time(2)-time(1));
    plot(diff_time*1e3,diff_sub_raw)
    hold on
    plot(diff_time*1e3,sub_dat_grad_smooth)
    plot(diff_time*1e3,pos_slope_mask*2e3)
    hold off
    ylim([min(diff_sub_raw),max(diff_sub_raw)])
    xlabel('time (ms)')
    ylabel('grad (v/s)')


    subplot(2,2,3)
    plot(sweep{ii}.pzt_smooth,sweep{ii}.pd_full_raw,'k')
    xlabel('pzt(v)')
    ylabel('pd (v)')
    hold on
    plot(sweep{ii}.pks.full.pzt,sweep{ii}.pks.full.pd,'xk','markersize',20)
    plot(sweep{ii}.pzt_smooth,sweep{ii}.pd_cmp_raw*cmp_multiplier_disp,'r')
    plot(sweep{ii}.pks.cmp.pzt,sweep{ii}.pks.cmp.pd*cmp_multiplier_disp,'rx','markersize',20);
    hold off
    pause(1e-5)
end


sweep{ii}.pks.all.pzt=[sweep{ii}.pks.cmp.pzt,sweep{ii}.pks.full.pzt];
sweep{ii}.pks.all.pd=[sweep{ii}.pks.cmp.pd,sweep{ii}.pks.full.pd];

sweep{ii}.pks.all.min_pzt_v_diff=min_diff(sweep{ii}.pks.all.pzt);

if sweep{ii}.pks.all.min_pzt_v_diff<peak_dist_min_pass
    fprintf('\nLASER IS NOT SINGLE MODE!!!\n')
    single_mode_vec(ii)=false;
     test=false;%give up if anything looks bad
else
    %fprintf('laser looks ok\n');
    single_mode_vec(ii)=true;
end

end
if ii>=iimax
    test=false;
end
end
fprintf('\n')
