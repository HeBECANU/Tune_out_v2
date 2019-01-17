
function ai_log_single_out=ai_log_single(args_single)
% single mode interogated over the entire range
%in
%  args_single.dir
%  args_single.fname
%  args_single.sfp.num_checks
%  args_single.sfp.thresh_cmp_peak
%  args_single.sfp.peak_dist_min_pass
%  args_single.pd.time_start
%  args_single.pd.time_stop

% outs
% ai_log_single_out.pd.mean
% ai_log_single_out.pd.std
% ai_log_single_out.pd.median
% ai_log_out.ok.sfp



%%load the data
path=strcat(args_single.dir,args_single.fname);
fid = fopen(path,'r');
raw_line = fread(fid,Inf,'*char')'; % this is the fastest string read method 3.5s for 20 files
fclose(fid);
ai_dat=jsondecode(raw_line);
samples= size(ai_dat.Data,2);
sr=ai_dat.sample_rate;
aquire_time=samples/sr;


sampl_start=max(1,ceil(args_single.pd.time_start*sr));
sampl_stop=min(samples,ceil(args_single.pd.time_stop*sr));
probe_pd_during_meas=ai_dat.Data(1,sampl_start:sampl_stop);

ai_log_single_out.pd.mean=mean(probe_pd_during_meas);
ai_log_single_out.pd.std=std(probe_pd_during_meas);
ai_log_single_out.pd.median=median(probe_pd_during_meas);


%% plot the anlog values
%as this function is not passed what the pd value should be it cant decide if it is failed for the
%args_single.plot.failed option
if args_single.plot.all || (args_single.plot.failed)      
    sfigure(1);
    set(gcf,'color','w')
    subplot(2,2,4)
    cla;
    plot((1:numel(probe_pd_during_meas))/sr,probe_pd_during_meas,'b')
    %yl=ylim;
    %xl=xlim;
%     hold on
%     line([xl(1),xl(2)],[args_single.pd.set,args_single.pd.set],'Color','k','LineWidth',3)
%     hold off
    %ylim([yl(1),yl(2)])
    ylabel('probe voltage')
    xlabel('time (s)')
    pause(1e-6)

end

%%
%%If the check if the laser is single mode
test_times=linspace(0,min(aquire_time-args_single.window_time,args_single.pd.time_stop),args_single.sfp.num_checks); 
pzt_division=0.2498; %set by the voltage divider box
peak_distance_thresh_cmp_full_min=0.3;%volts difference full-cmp peak in order to be considered new peak

test_sm_while=true; %intialize while loop flag, give up with one bad detection
jj=0;
jjmax=numel(test_times);
sweep=[];
single_mode_vec=false(1,jjmax);
sfp_pzt_raw=ai_dat.Data(3,:);
sfp_pd_raw=ai_dat.Data(2,:);
jj=0;
while test_sm_while
    jj=jj+1;
    time_start=test_times(jj);
    sampl_start=max(1,1+floor(time_start*sr));
    sampl_stop=min(samples,sampl_start+ceil(args_single.window_time*sr));
    if sampl_stop-sampl_start>=args_single.window_time*sr %check that we have enough points to work with
        sub_ptz_raw=sfp_pzt_raw(sampl_start:sampl_stop)/pzt_division;
        kernel=gausswin(ceil(4*args_single.pzt_volt_smothing_time*sr),4);
        kernel=kernel/sum(kernel(:));%normalize
        sub_dat_smooth=conv(sub_ptz_raw,kernel,'same');
        sub_dat_grad_smooth=diff(sub_dat_smooth)*sr;
        pos_slope_mask=sub_dat_grad_smooth>0;
        %now we want to select the first positive going sgement that is long enough
        x=pos_slope_mask';
        % Find the transitions
        down= strfind(x',[1 0]);
        ups= strfind(x',[0 1]);
        %we take the ramp to investigate as the between the first rising edge of the mask and the next
        %falling edge
        sweep{jj}.start_idx=ups(1);
        %find the next falling edge
        sweep{jj}.stop_idx=down(find(down>sweep{jj}.start_idx,1)); %fuck so elegant, fuck this is so much beter than LabView
        sweep{jj}.pzt_raw=sub_ptz_raw(sweep{jj}.start_idx:sweep{jj}.stop_idx);
        sweep{jj}.pzt_smooth=sub_dat_smooth(sweep{jj}.start_idx:sweep{jj}.stop_idx);
        sweep{jj}.pd_full_raw=sfp_pd_raw(sampl_start+sweep{jj}.start_idx:sampl_start+sweep{jj}.stop_idx);
        sweep{jj}.pd_cmp_raw=ai_dat.Data(4,sampl_start+sweep{jj}.start_idx:sampl_start+sweep{jj}.stop_idx);
        sweep{jj}.pd_cmp_raw=sweep{jj}.pd_cmp_raw-median(sweep{jj}.pd_cmp_raw);
        %threshold out all data that is 7% the max value for the peak finding thresh
        thesh=0.07*max(sweep{jj}.pd_full_raw);
        [pks_val,locs] = findpeaks(sweep{jj}.pd_full_raw,'MinPeakHeight',thesh);
        sweep{jj}.pks.full.pd=pks_val;
        sweep{jj}.pks.full.pzt=sweep{jj}.pzt_smooth(locs);
        [pks_val,locs] = findpeaks(sweep{jj}.pd_cmp_raw,'MinPeakHeight',args_single.sfp.thresh_cmp_peak);
        sweep{jj}.pks.cmp.pd=pks_val;
        sweep{jj}.pks.cmp.pzt=sweep{jj}.pzt_smooth(locs);
        %find the distance in pzt voltage to the nearest full peak
        min_dist_cmp_to_full=arrayfun(@(x) min(abs(x-sweep{jj}.pks.full.pzt)),sweep{jj}.pks.cmp.pzt);
        %peak_distance_thresh_cmp_full_min
        sweep{jj}.pks.cmp.pzt=sweep{jj}.pks.cmp.pzt(min_dist_cmp_to_full>peak_distance_thresh_cmp_full_min);
        sweep{jj}.pks.cmp.pd=sweep{jj}.pks.cmp.pd(min_dist_cmp_to_full>peak_distance_thresh_cmp_full_min);
        sweep{jj}.pks.all.pzt=[sweep{jj}.pks.cmp.pzt,sweep{jj}.pks.full.pzt];
        sweep{jj}.pks.all.pd=[sweep{jj}.pks.cmp.pd,sweep{jj}.pks.full.pd];
        sweep{jj}.pks.all.min_pzt_v_diff=min_diff(sweep{jj}.pks.all.pzt);

        if sweep{jj}.pks.all.min_pzt_v_diff<args_single.sfp.peak_dist_min_pass
            fprintf('\nLASER IS NOT SINGLE MODE!!!\n%04i',0')
            single_mode_vec(jj)=false;%give up if anything looks bad
        else
           single_mode_vec(jj)=true;

        end

        if args_single.plot.all || (args_single.plot.failed && ~single_mode_vec(jj))
            title_strs={'Failed','OK'};
            sfigure(1);
            subplot(2,2,1)
            cla;
             set(gcf,'color','w')
            title('smoothing pzt v')
            time=(1:numel(sub_ptz_raw))/sr;
            plot(time*1e3,sub_ptz_raw,'r')
            hold on
            plot(time*1e3,sub_dat_smooth,'k')
            hold off
            xlabel('time (ms)')
            ylabel('volts (v)')

            subplot(2,2,2)
            cla;
            %select a single scan
            %to do this we will mask out the positive slope
            diff_sub_raw=diff(sub_ptz_raw)*sr;
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
            cla;
            plot(sweep{jj}.pzt_smooth,sweep{jj}.pd_full_raw,'k')
            xlabel('pzt(v)')
            ylabel('pd (v)')
            hold on
            plot(sweep{jj}.pks.full.pzt,sweep{jj}.pks.full.pd,'xk','markersize',20)
            plot(sweep{jj}.pzt_smooth,sweep{jj}.pd_cmp_raw*args_single.cmp_multiplier_disp,'r')
            plot(sweep{jj}.pks.cmp.pzt,sweep{jj}.pks.cmp.pd*args_single.cmp_multiplier_disp,'rx','markersize',20);
            hold off
            title(title_strs(single_mode_vec(jj)+1))

            pause(1e-6)
        end %end plots
    end %enough points to work with

if jj>=jjmax || ~single_mode_vec(jj) %if something is not single mode do not continue
    test_sm_while=false;
end 
end%end while loop over checking mode
if sum(~single_mode_vec)==0
    ai_log_single_out.single_mode=true;
else
    ai_log_single_out.single_mode=false;
end

    
    
end