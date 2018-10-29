function ai_log_out=ai_log_import(anal_opts,data)
%a simple wrapper for the below ai_log_import that uses the matlab function cache
cache_opts=[];
cache_opts.verbose=3;
cache_opts.force_cache_load=anal_opts.force_load_save;
anal_opts=rmfield(anal_opts,'force_load_save');
cache_opts.force_recalc=anal_opts.force_reimport;
if ~cache_opts.force_recalc
    anal_opts=rmfield(anal_opts,'force_reimport');
end

%limit the scope but retain the structure
data_sub=[];
data_sub.mcp_tdc.shot_num=data.mcp_tdc.shot_num;
data_sub.mcp_tdc.time_create_write=data.mcp_tdc.time_create_write;

outputs=function_cache(cache_opts,@ai_log_import_core,{anal_opts,data_sub});
ai_log_out=outputs{1};
end

function ai_log_out=ai_log_import_core(anal_opts,data)
%ai_log_import - imports analog input log and checks if the probe beam pd signal is ok and that the laser is single
%mode by measuring the distance (in pzt voltage) between the scanning FP pd peaks
%results are placed into a convenient strucure with mat file cache
%will load data from a cashed version if anal_opts has not changed
%the output is a well aranged structure to be added into the data structure
%the data structure is not included in the save to prevent double saving of the data
%at the end of the import or load cashe the approapriate feilds are added to data
%Only checks if sm if pd is ok/all previous checks of sm are ok
%
% Syntax:  [data,import_opts]=import_data(import_opts)
%

%inputs

%output
    %ai_log_out struct
    %.file
% Known BUGS/ Possible Improvements
%   -need more dam speed!
%       -reading and jsondecode are the main problems the rest is very fast
%   -Improved documentation
%   -



% Author: Bryce Henson
% email: Bryce.Henson[a circle]live.com  %YOU MUST INCLUDE '[matlab][ai_log_import]' in the subject line OR I WILL NOT REPLY
% Last revision:2018-09-30

%------------- BEGIN USER VAR --------------
%estimate of the sfp scan time,used to set the window and the smoothing
args_single=[];
args_single.cmp_multiplier_disp=50; %multiplier to display the compressed data better
args_single.window_time=anal_opts.scan_time*2.1;
args_single.pzt_volt_smothing_time=anal_opts.scan_time/100;
%------------- END USER VAR --------------

%------------- BEGIN CODE --------------

dir_read=dir([anal_opts.dir,anal_opts.log_name,'*.txt']);
ai_log_out.file_names={dir_read.name};

%now find the shot that this corresponds to
%time of start shot on TDC_comp
%use write so files are portable
time_start_tdc_comp=data.mcp_tdc.time_create_write(:,2)-anal_opts.trig_dld-anal_opts.dld_aquire;

%the number of shots in the mcp_tdc struct, note not ness the same as the number of ai_logs
shots_tdc=size(data.mcp_tdc.shot_num,2);
%GUILTY UNTILL PROVEN INOCENT !!!!!!!!!!!!!!!!!
ai_log_out.ok.reg_pd=false(shots_tdc,1);  %pd regulaton
ai_log_out.ok.sfp=false(shots_tdc,1);     %scanning FP check
dld_files=numel(data.mcp_tdc.shot_num);

%set up for the ai_log_single 
args_single.time_start_tdc_comp=time_start_tdc_comp;
args_single.log_name=anal_opts.log_name;
args_single.dir=anal_opts.dir;
args_single.trig_ai_in=anal_opts.trig_ai_in;
args_single.time_match_valid=anal_opts.time_match_valid;
args_single.mdp_shot_num=data.mcp_tdc.shot_num ;
args_single.pd=anal_opts.pd;
args_single.plot=anal_opts.plot;
args_single.sfp=anal_opts.sfp;

cache_opts=[];
cache_opts.verbose=0;
cache_opts.mock_working_dir=anal_opts.dir;
cache_opts.path_directions={1,'dir'};
cache_opts.clean_cache=false;
if isfield(anal_opts,'force_reimport')
    %cache_opts.force_recalc=true;
end


iimax=size(ai_log_out.file_names,2); %the number of ai logs that have been identified
%initalize outputs
ai_log_out.ok.reg_pd=false(dld_files,1);
ai_log_out.ok.sfp=false(dld_files,1);
%loop over all the ai_logs
fprintf('processing ai log for files %04u:%04u',iimax,0)
for ii=1:iimax
    fprintf('\b\b\b\b%04i',ii)
    args_single.fname=ai_log_out.file_names{ii};
    cout=function_cache(cache_opts,@ai_log_single,{args_single});
    ai_log_single_out=cout{1};
    ai_log_out.ok.reg_pd(ai_log_single_out.shot_idx)=ai_log_single_out.reg_pd;
    ai_log_out.ok.sfp(ai_log_single_out.shot_idx)=ai_log_single_out.single_mode;
    
end %loop over files
fprintf('Done\n')

end%function


function ai_log_single_out=ai_log_single(args_single)
% outs
% ai_log_single_out.single_mode
% ai_log_single_out.reg_pd
%ai_log_single_out.shot_idx
%ai_log_single_out.shot_num

fname=args_single.fname;
time_iso_str=erase(erase(fname,'log_analog_in_'),'.txt');
time_start_tdc_comp=args_single.time_start_tdc_comp;

path=strcat(args_single.dir,fname);
%bit of a hack to get data_tcreate to work which was set up to take a path+num format \d123.txt
time_posix_ai_log_create_write=data_tcreate([args_single.dir,args_single.log_name],time_iso_str);
fid = fopen(path,'r');
raw_line = fread(fid,Inf,'*char')'; % this is the fastest string read method 3.5s for 20 files
fclose(fid);

time_posix_fname=posixtime(datetime(time_iso_str,'InputFormat','yyyyMMdd''T''HHmmss'));
ai_dat=jsondecode(raw_line);
samples= size(ai_dat.Data,2);
sr=ai_dat.sample_rate;
aquire_time=samples/sr;
%initalize the logic as failed
ai_log_single_out.reg_pd=false;
ai_log_single_out.single_mode=false;
ai_log_single_out.shot_idx=nan;

%time of the start on the BEC comp
%for this to work the clock sync need to be decent
time_start_bec_comp=time_posix_ai_log_create_write(2)-args_single.trig_ai_in-aquire_time;
[time_nearest_tdc_start,idx_nearest_shot]=closest_value(time_start_tdc_comp,time_start_bec_comp);
%should not process if not near a shot
if abs(time_nearest_tdc_start-time_start_bec_comp)>args_single.time_match_valid
     fprintf(2,'\nnearest tdc file is too far away at %.2f s\n%04u',time_nearest_tdc_start-time_start_bec_comp,0)
else
    ai_log_single_out.shot_idx=idx_nearest_shot; %index in the mcp_tdc arrays of this shot

    ai_log_single_out.shot_num=args_single.mdp_shot_num(idx_nearest_shot); %number of shot eg d54.txt
    ai_log_single_out.tdiff=time_nearest_tdc_start-time_start_bec_comp;
    %output all the times mainly as a diagnostic
    ai_log_single_out.times.create_ai_log=time_posix_ai_log_create_write(1); 
    ai_log_single_out.times.fname_ai_log=time_posix_fname;

    sampl_start=max(1,ceil(args_single.pd.time_start*sr));
    sampl_stop=min(samples,ceil(args_single.pd.time_stop*sr));
    probe_pd_during_meas=ai_dat.Data(1,sampl_start:sampl_stop);

    ai_log_single_out.pd.mean=mean(probe_pd_during_meas);
    ai_log_single_out.pd.std=std(probe_pd_during_meas);

    if abs(ai_log_single_out.pd.mean-args_single.pd.set)>args_single.pd.diff_thresh
        fprintf('\nprobe beam pd value wrong!!!!!!!!!\n',0)
        fprintf('avg val %.2f set value %.2f \n04u%',ai_log_single_out.pd.mean,args_single.pd.set,0)
    elseif ai_log_single_out.pd.std>args_single.pd.std_thresh
        fprintf('\nprobe beam pd noisy!!!!!!!!!\n%04u',0)
        fprintf('std %.2f thresh value %.2f \n04u%',ai_log_single_out.pd.std,args_single.pd.std_thresh,0)
    else
        ai_log_single_out.reg_pd=true;
    end

    if args_single.plot.all || (args_single.plot.failed && ~ai_log_single_out.reg_pd)      
        sfigure(1);
        subplot(2,2,4)
        title_strs={'Failed','OK'};
        set(gcf,'color','w')
        plot((1:numel(probe_pd_during_meas))/sr,probe_pd_during_meas,'b')
        yl=ylim;
        xl=xlim;
        hold on
        line([xl(1),xl(2)],[args_single.pd.set,args_single.pd.set],'Color','k','LineWidth',3)
        hold off
        ylim([yl(1),yl(2)])
        ylabel('probe voltage')
        xlabel('time (s)')
        title(title_strs(ai_log_single_out.reg_pd+1))
        pause(1e-6)

    end
end
%%
%%If the pd is ok then lets check if the laser is single mode
if ai_log_single_out.reg_pd
    test_times=linspace(0,min(aquire_time-args_single.window_time,args_single.pd.time_stop),args_single.sfp.num_checks); 

    anal_opts.sfp.thresh_cmp_peak=20e-3;%volts
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
                fprintf('\nLASER IS NOT SINGLE MODE!!!\n%03u',0')
                single_mode_vec(jj)=false;%give up if anything looks bad
            else
               single_mode_vec(jj)=true;

            end

            if args_single.plot.all || (args_single.plot.failed && ~single_mode_vec(jj))
                title_strs={'Failed','OK'};
                sfigure(1);
                subplot(2,2,1)
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
                plot(sweep{jj}.pzt_smooth,sweep{jj}.pd_full_raw,'k')
                xlabel('pzt(v)')
                ylabel('pd (v)')
                hold on
                plot(sweep{jj}.pks.full.pzt,sweep{jj}.pks.full.pd,'xk','markersize',20)
                plot(sweep{jj}.pzt_smooth,sweep{jj}.pd_cmp_raw*cmp_multiplier_disp,'r')
                plot(sweep{jj}.pks.cmp.pzt,sweep{jj}.pks.cmp.pd*cmp_multiplier_disp,'rx','markersize',20);
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
end%if pd test passed
    
    
end