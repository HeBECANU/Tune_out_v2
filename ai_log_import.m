function [data,ai_log_out]=ai_log_import(anal_opts,data)
%ai_log_import - imports analog input log and checks if the probe beam pd signal is ok and that the laser is single
%mode by measuring the distance (in pzt voltage) between the scanning FP pd peaks
%results are placed into a convineint strucure with mat file cache
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

%------------- BEGIN CODE --------------
%begin user var
scan_time=25e-3; %non citital, just needs to be long enough to capture a few scans of the FP
%end user var




import_logs=true;
%should improve this and save it in the data directory
if  isfile('import_ai_log_save.mat') && (anal_opts.force_load_save ||  ~anal_opts.force_reimport )
    load('import_ai_log_save.mat','anal_opts_old')
    anal_opts_old.force_reimport=anal_opts.force_reimport; %prevents reimport after import_opts.force_reimport changed 1->0
    if isequal(anal_opts_old,anal_opts) || anal_opts.force_load_save
        import_logs=false; 
        fprintf('anal_opts the same loading old data...')
        load('import_ai_log_save.mat','ai_log_out','anal_opts_old')
        fprintf('Done\n')

    end
end




if import_logs
    dir_read=dir([anal_opts.dir,anal_opts.log_name,'*.txt']);
    ai_log_out.ai_log.file_names={dir_read.name};
    
    %now find the shot that this corresponds to
    %time of start shot on TDC_comp
    time_start_tdc_comp=data.mcp_tdc.time_create_write(:,1)-anal_opts.trig_dld;

    %the number of shots in the mcp_tdc struct, note not ness the same as the number of ai_logs
    shots_tdc=size(data.mcp_tdc.shot_num,2);
    %GUILTY UNTILL PROVEN INOCENT !!!!!!!!!!!!!!!!!
    ai_log_out.mcp_tdc.probe.ok.reg_pd=false(shots_tdc,1);  %pd regulaton
    ai_log_out.mcp_tdc.probe.ok.sfp=false(shots_tdc,1);     %scanning FP check
    
    iimax=size(ai_log_out.ai_log.file_names,2); %the number of ai logs that I have identified
    %loop over all the ai_logs
    fprintf('processing ai log for files %04u:%04u',iimax,0)
    for ii=1:iimax
        fprintf('\b\b\b\b%04i',ii)
        fname=ai_log_out.ai_log.file_names{ii};
        time_iso_str=erase(erase(fname,'log_analog_in_'),'.txt');
        path=strcat(anal_opts.dir,fname);
        %bit of a hack to get data_tcreate to work which was set up to take a path+num format \d123.txt
        time_posix_ai_log_create_write=data_tcreate([anal_opts.dir,anal_opts.log_name],time_iso_str);
        fid = fopen(path,'r');
        %reading is one of the slowest parts of the code
        %ai_log_file_cells=textscan(fid,'%s','Delimiter','\n'); %5.9s for 20 files
        %line=fgets(fid); %11.7s
        line = fread(fid,Inf,'*char')'; %3.5s for 20 files
        fclose(fid);

        ai_dat=jsondecode(line);
        samples= size(ai_dat.Data,2);
        sr=ai_dat.sample_rate;
        aquire_time=samples/sr;

        time_posix_fname=posixtime(datetime(time_iso_str,'InputFormat','yyyyMMdd''T''HHmmss'));
        %time diff between the file name and the properties created time
        %this is about the cycle duration as the fname is created after the last file is written and it is setting up
        %input
%         time_diff_name_prop=time_posix_fname-time_posix_ai_log_create_write(1);
%         if time_diff_name_prop>2
%             warning('time difference between file name and file creation is very large\n')
%         end


        %time of the start on the BEC comp
        %for this to work the clock sync need to be decent
        time_start_bec_comp=time_posix_ai_log_create_write(1)-anal_opts.trig_ai_in-aquire_time;
        [time_nearest_tdc_start,idx_nearest_shot]=closest_value(time_start_tdc_comp,time_start_bec_comp);
        %should not process if not near a shot
        ai_log_out.ai_log.shot_idx(ii)=idx_nearest_shot; %index in the mcp_tdc arrays of this shot
        ai_log_out.ai_log.shot_num(ii)=data.mcp_tdc.shot_num(idx_nearest_shot); %number of shot eg d54.txt
        ai_log_out.ai_log.tdiff(ii)=time_nearest_tdc_start-time_start_bec_comp;
        %output all the times mainly as a diagnostic
        ai_log_out.ai_log.times.create_ai_log=time_posix_ai_log_create_write(1); 
        ai_log_out.ai_log.times.fname_ai_log=time_posix_fname;

        sampl_start=max(1,ceil(anal_opts.pd.time_start*sr));
        sampl_stop=min(samples,ceil(anal_opts.pd.time_stop*sr));
        probe_pd_during_meas=ai_dat.Data(1,sampl_start:sampl_stop);

        ai_log_out.ai_log.pd.mean(ii)=mean(probe_pd_during_meas);
        ai_log_out.ai_log.pd.std(ii)=std(probe_pd_during_meas);

        if abs(ai_log_out.ai_log.pd.mean(ii)-anal_opts.pd.set)>anal_opts.pd.diff_thresh
            fprintf('\nprobe beam pd value wrong!!!!!!!!!\n%04u',0)
        elseif ai_log_out.ai_log.pd.std(ii)>anal_opts.pd.std_thresh
            fprintf('\nprobe beam pd noisy!!!!!!!!!\n%04u',0)
        else
            ai_log_out.mcp_tdc.probe.ok.reg_pd(ii)=true;
        end
        sfigure(2);
        clf
        plot(probe_pd_during_meas)





        %%
        %%If the pd is ok then lets check if the laser is sm
        if ai_log_out.mcp_tdc.probe.ok.reg_pd(ii)
            test_times=linspace(0,min(aquire_time-scan_time,anal_opts.pd.time_stop),anal_opts.sfp.num_checks); 

            anal_opts.sfp.thresh_cmp_peak=20e-3;%volts
            pzt_division=0.2498; %set by the voltage divider box
            peak_distance_thresh_cmp_full_min=0.3;%volts difference full-cmp peak in order to be considered new peak


            anal_opts.sfp.plot_all=false;
            anal_opts.sfp.plot_failed=true;

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
                sampl_start=max(1,ceil(time_start*sr));
                sampl_stop=min(samples,sampl_start+ceil(scan_time*sr));
                if sampl_stop-sampl_start>=scan_time*sr %check that we have enough points to work with
                    sub_dat_raw=sfp_pzt_raw(sampl_start:sampl_stop)/pzt_division;
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

                    sweep{jj}.start_idx=ups(1);
                    %find the next falling edge
                    sweep{jj}.stop_idx=down(find(down>sweep{jj}.start_idx,1)); %fuck so elegant, fuck this is so much beter than LabView

                    sweep{jj}.pzt_raw=sub_dat_raw(sweep{jj}.start_idx:sweep{jj}.stop_idx);
                    sweep{jj}.pzt_smooth=sub_dat_smooth(sweep{jj}.start_idx:sweep{jj}.stop_idx);
                    sweep{jj}.pd_full_raw=sfp_pd_raw(sampl_start+sweep{jj}.start_idx:sampl_start+sweep{jj}.stop_idx);
                    sweep{jj}.pd_cmp_raw=ai_dat.Data(4,sampl_start+sweep{jj}.start_idx:sampl_start+sweep{jj}.stop_idx);
                    sweep{jj}.pd_cmp_raw=sweep{jj}.pd_cmp_raw-median(sweep{jj}.pd_cmp_raw);
                    %threshold out all data that is 7% the max value for the peak finding thresh
                    thesh=0.07*max(sweep{jj}.pd_full_raw);

                    cmp_multiplier_disp=50; %multiplier to display the compressed data better


                    [pks_val,locs] = findpeaks(sweep{jj}.pd_full_raw,'MinPeakHeight',thesh);
                    sweep{jj}.pks.full.pd=pks_val;
                    sweep{jj}.pks.full.pzt=sweep{jj}.pzt_smooth(locs);

                    [pks_val,locs] = findpeaks(sweep{jj}.pd_cmp_raw,'MinPeakHeight',anal_opts.sfp.thresh_cmp_peak);
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

                    if sweep{jj}.pks.all.min_pzt_v_diff<anal_opts.sfp.peak_dist_min_pass
                        fprintf('\nLASER IS NOT SINGLE MODE!!!\n%03u',0')
                        single_mode_vec(jj)=false;
                         test_sm_while=false;%give up if anything looks bad
                    else
                        single_mode_vec(jj)=true;
                        %fprintf('laser looks ok\n');

                    end


                    if anal_opts.sfp.plot_all || (anal_opts.sfp.plot_failed && ~single_mode_vec(jj))
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
                        plot(sweep{jj}.pzt_smooth,sweep{jj}.pd_full_raw,'k')
                        xlabel('pzt(v)')
                        ylabel('pd (v)')
                        hold on
                        plot(sweep{jj}.pks.full.pzt,sweep{jj}.pks.full.pd,'xk','markersize',20)
                        plot(sweep{jj}.pzt_smooth,sweep{jj}.pd_cmp_raw*cmp_multiplier_disp,'r')
                        plot(sweep{jj}.pks.cmp.pzt,sweep{jj}.pks.cmp.pd*cmp_multiplier_disp,'rx','markersize',20);
                        hold off
                        pause(1e-5)
                    end %end plots
                end %enough points to work with
                
            if jj>=jjmax || ~single_mode_vec(jj) %if something is not single mode do not continue
                test_sm_while=false;
            end 
            end%end while loop over checking mode
            if sum(~single_mode_vec)==0
                ai_log_out.mcp_tdc.probe.ok.sfp(ii)=true;
            end
        end%if pd test passed

    end %loop over files
    fprintf('\b\b\b\b...Done\nSaving Data...')
    anal_opts_old=anal_opts;
    save('import_ai_log_save.mat','ai_log_out','anal_opts_old','-v7.3') %,'-nocompression'
    fprintf('Done\n')
end%import data
end%function