clear all
close all

% TODO: Correct for PD offset
% TODO: Check zero-peak files
% TODO: Using calibration offset, compute mean power ratio per file
% Then compile and interpret (i.e. physics)

treshold = 1.5;
zero_offset = 0;
gamma_multiplier = 20;
peak_width = 1.5e-4;
logscale = false;
lorentz_fit = false;
plot_all = false;

args.dir = 'C:\Data\blue_sfp\red_on\';

anal_opts.log_name='log_analog_in_';

anal_opts.scan_time=14e-3;  %estimate of the sfp scan time,used to set the window and the smoothin
args.pzt_volt_smothing_time=anal_opts.scan_time/100;

dir_read=dir([args.dir,anal_opts.log_name,'*.txt']);
file_names={dir_read.name};
num_files = numel(file_names);


args.fname = 'log_analog_in_20181224T143228.txt';

profile on
dir_stats = zeros(num_files,2);
for pp = 1:num_files
    pp
    args.fname = file_names{pp};

    %%load the data
    path=strcat(args.dir,args.fname);
    fid = fopen(path,'r');
    raw_line = fread(fid,Inf,'*char')'; % this is the fastest string read method 3.5s for 20 files
    fclose(fid);
    ai_dat=jsondecode(raw_line);
    samples= size(ai_dat.Data,2);
    sr=ai_dat.sample_rate;
    acquire_time=samples/sr;
    T = linspace(0,acquire_time, samples);

    %ai_dat format:
    %   1: probe beam PD; low, because attenuated by SFP?
    % sfp_pzt_raw=ai_dat.Data(3,:);
    % sfp_pd_raw=ai_dat.Data(2,:);
    % sweep{jj}.pd_cmp_raw=ai_dat.Data(4,sampl_start+sweep{jj}.start_idx:sampl_start+sweep{jj}.stop_idx);


    sampl_start = 1;
    demo_max = 4e3;

    pzt_division=0.2498;
    sfp_pzt_raw=ai_dat.Data(3,:);
    sub_ptz_raw=sfp_pzt_raw(sampl_start:demo_max)/pzt_division;
    kernel=gausswin(ceil(4*args.pzt_volt_smothing_time*sr),4);
    kernel=kernel/sum(kernel(:));%normalize
    sub_dat_smooth=conv(sub_ptz_raw,kernel,'same');
    sub_dat_grad_smooth=diff(sub_dat_smooth)*sr;
    pos_slope_mask=sub_dat_grad_smooth>0;
    pos_scan_data = ai_dat.Data(2,pos_slope_mask);

    %now we want to select the first positive going sgement that is long enough
    x=pos_slope_mask';
    % Find the transitions
    downs = strfind(x',[1 0]);
    ups = strfind(x',[0 1]);


    pd_raw = ai_dat.Data(2,:);

    %Use only complete scans
    T_min = min(ups);
    T_max = max(downs);
    downs = downs(downs > T_min);
    ups = ups(ups < T_max);
    borders = [ups',downs'];
    num_scans = size(borders,1);



    figure();
    clf;
    
    file_ratios = zeros(num_scans,2);
    
    for ii=1:num_scans

        idx_lims = borders(ii,:);
        sweep = pd_raw(idx_lims(1):idx_lims(2));
        T_sweep = T(idx_lims(1):idx_lims(2));
        [pks_val,locs] = findpeaks(sweep,'MinPeakHeight',treshold);


        if ~isempty(locs)

            hwin_size = floor(0.5*median(diff(locs)));
            num_peaks = length(locs);
            scan_ratios = zeros(num_peaks, 1);
            
            for jj = 1:num_peaks
                p_cent = locs(jj);

                idx_lims = [(max(1,p_cent-hwin_size)), min(length(sweep),p_cent + hwin_size)];
                T_win = T_sweep(idx_lims(1):idx_lims(2));
                p_edge = T_sweep(idx_lims);

                if diff(idx_lims) > 1.5*hwin_size % ensure all or most of a scan is taken
        %             calc_pow_ratio()
                    scan = sweep(idx_lims(1):idx_lims(2));
                    scan = scan - zero_offset;
                    T_cen = T_sweep(p_cent);

                    if lorentz_fit
                        beta0 = [pks_val(jj),1e-5,T_cen,zero_offset];
                        mdl = fitnlm(T_win,scan,@Lorentzian,beta0);
                        coefs = mdl.Coefficients.Estimate;
                        gamma = coefs(2);
                        T_cen = coefs(3);
    %                     line_lim = gamma_multiplier*gamma;
                    else
    %                     line_lim = peak_width;
                    end
                    
                    
                    line_lim = peak_width;
                    line_lims = [T_cen - line_lim,T_cen + line_lim];
                    peak_mask = T_win > line_lims(1) & T_win< line_lims(2);
                    dT = mean(diff(T_win)); %Not correct but need a quick fix
                    peak_power = sum(scan(peak_mask)*dT);
                    back_power = sum(scan(~peak_mask)*dT);
                    T_fit = linspace(min(T_win),max(T_win),2000);
                    scan_ratios(jj) = back_power/peak_power;
                    
                    if plot_all
                        subplot(2,1,1)
                        plot(T_win,scan,'k.')
                        hold on
                        if logscale
                            plot(T_fit,Lorentzian(coefs,T_fit),'b')
                            set(gca,'Yscale','log')
                        end
                        pl =  gca;
                        ylim = pl.YLim;
                        plot([p_edge;p_edge],[ylim(1),ylim(1);ylim(2),ylim(2)],'r-.')
                        plot([line_lims;line_lims],[ylim(1),ylim(1);ylim(2),ylim(2)],'g-.')
                        xlim([0,max(T_sweep)])

                        subplot(4,1,3)
                        plot(T_cen,peak_power,'kx')
                        hold on

                        subplot(4,1,4)
                        plot(T_cen,scan_ratios(jj),'kx')
                        hold on

                        xlim([0,max(T_sweep)])
                    end

                end


            end %loop over peaks

        end %isempty(locs)

        
    end %loop over scans
    
    if plot_all
        subplot(2,1,1)
        title('SFP scan')
        xlabel('time [s]')
        ylabel('PD voltage')
        legend('Lorentzian fit','Data','Peak borders')

        subplot(4,1,3)
        xlim([0,max(T_sweep)])
        title('Peak power')
        xlabel('Time')
        ylabel('Power')

        subplot(4,1,4)
        title('Background/peak power')
        xlabel('time [s]')
        ylabel('Power ratio')
        suptitle(sprintf('%f', dir_stats(pp,1)));
    end
    dir_stats(pp,:) = [mean(scan_ratios),std(scan_ratios)];
    

end
dir_stats
profile off
profile viewer
% Breaking up by FSR

    function y = Lorentzian(p,x)
        gamma = p(2);
        amp = 2*p(1)*gamma;
        x0 = p(3);
        dc = p(4);
        y = amp*(1/pi)*(0.5*gamma)./((x-x0).^2 + (0.5*gamma).^2)+p(4);
    end


%% Bryce scan detection method:

% args.trig_ai_in=anal_opts.trig_ai_in;
% args.plot=anal_opts.plot;
% args.sfp=anal_opts.sfp;
% args.pd.time_start=anal_opts.pd.time_start;
% args.pd.time_stop=anal_opts.pd.time_stop;
% % args.plot.all = true;
% args.pd.time_start = 0;
% args.pd.time_stop = 4;
% args.window_time = 2;
% 
% 
% anal_opts.trig_ai_in=20;
% anal_opts.force_load_save=false;
% anal_opts.log_name='log_analog_in_';
% anal_opts.acquire_time=3;
% anal_opts.pd.diff_thresh=0.1;
% anal_opts.pd.std_thresh=0.1;
% anal_opts.pd.time_start=0.2;
% anal_opts.pd.time_stop=2;
% anal_opts.sfp.num_checks=10; %how many places to check that the laser is single mode
% anal_opts.sfp.thresh_cmp_peak=5e-2; %theshold on the compressed signal to be considered a peak
% anal_opts.sfp.peak_dist_min_pass=4.5;%minimum (min difference)between peaks for the laser to be considered single mode
% anal_opts.plot.all=false;
% anal_opts.plot.failed=false;
% anal_opts.time_match_valid=5; %how close the predicted start of the shot is to the actual


% args.fname = 'log_analog_in_20181224T143228.txt';


% sampl_start=max(1,ceil(args.pd.time_start*sr));
% sampl_stop=min(samples,ceil(args.pd.time_stop*sr));
% probe_pd_during_meas=ai_dat.Data(1,sampl_start:sampl_stop);
% 
% ai_log_single_out.pd.mean=mean(probe_pd_during_meas);
% ai_log_single_out.pd.std=std(probe_pd_during_meas);
% ai_log_single_out.pd.median=median(probe_pd_during_meas);



%%
%%If the check if the laser is single mode
% test_times=linspace(0,min(acquire_time-args.window_time,args.pd.time_stop),args.sfp.num_checks); 
% pzt_division=0.2498; %set by the voltage divider box
% peak_distance_thresh_cmp_full_min=0.3;%volts difference full-cmp peak in order to be considered new peak
% 
% test_sm_while=true; %intialize while loop flag, give up with one bad detection
% jjmax=numel(test_times);
% sweep=[];
% single_mode_vec=false(1,jjmax);
% sfp_pzt_raw=ai_dat.Data(3,:);
% sfp_pd_raw=ai_dat.Data(2,:);
% jj=0;
% while test_sm_while
%     jj=jj+1;
%     time_start=test_times(jj);
%     sampl_start=max(1,1+floor(time_start*sr));
%     sampl_stop=min(samples,sampl_start+ceil(args.window_time*sr));
%     if sampl_stop-sampl_start>=args.window_time*sr %check that we have enough points to work with
%         sub_ptz_raw=sfp_pzt_raw(sampl_start:sampl_stop)/pzt_division;
%         kernel=gausswin(ceil(4*args.pzt_volt_smothing_time*sr),4);
%         kernel=kernel/sum(kernel(:));%normalize
%         sub_dat_smooth=conv(sub_ptz_raw,kernel,'same');
%         sub_dat_grad_smooth=diff(sub_dat_smooth)*sr;
%         pos_slope_mask=sub_dat_grad_smooth>0;
%         %now we want to select the first positive going sgement that is long enough
%         x=pos_slope_mask';
%         % Find the transitions
%         down= strfind(x',[1 0]);
%         ups= strfind(x',[0 1]);
%         %we take the ramp to investigate as the between the first rising edge of the mask and the next
%         %falling edge
%         sweep{jj}.start_idx=ups(1);
%         %find the next falling edge
%         sweep{jj}.stop_idx=down(find(down>sweep{jj}.start_idx,1)); %fuck so elegant, fuck this is so much beter than LabView
%         sweep{jj}.pzt_raw=sub_ptz_raw(sweep{jj}.start_idx:sweep{jj}.stop_idx);
%         sweep{jj}.pzt_smooth=sub_dat_smooth(sweep{jj}.start_idx:sweep{jj}.stop_idx);
%         sweep{jj}.pd_full_raw=sfp_pd_raw(sampl_start+sweep{jj}.start_idx:sampl_start+sweep{jj}.stop_idx);
%         sweep{jj}.pd_cmp_raw=ai_dat.Data(4,sampl_start+sweep{jj}.start_idx:sampl_start+sweep{jj}.stop_idx);
%         sweep{jj}.pd_cmp_raw=sweep{jj}.pd_cmp_raw-median(sweep{jj}.pd_cmp_raw);
%         %threshold out all data that is 7% the max value for the peak finding thresh
%         thesh=0.07*max(sweep{jj}.pd_full_raw);
%         [pks_val,locs] = findpeaks(sweep{jj}.pd_full_raw,'MinPeakHeight',thesh);
%         sweep{jj}.pks.full.pd=pks_val;
%         sweep{jj}.pks.full.pzt=sweep{jj}.pzt_smooth(locs);
%         [pks_val,locs] = findpeaks(sweep{jj}.pd_cmp_raw,'MinPeakHeight',args.sfp.thresh_cmp_peak);
%         sweep{jj}.pks.cmp.pd=pks_val;
%         sweep{jj}.pks.cmp.pzt=sweep{jj}.pzt_smooth(locs);
%         %find the distance in pzt voltage to the nearest full peak
% %         min_dist_cmp_to_full=arrayfun(@(x) min(abs(x-sweep{jj}.pks.full.pzt)),sweep{jj}.pks.cmp.pzt);
% %         %peak_distance_thresh_cmp_full_min
% %         sweep{jj}.pks.cmp.pzt=sweep{jj}.pks.cmp.pzt(min_dist_cmp_to_full>peak_distance_thresh_cmp_full_min);
% %         sweep{jj}.pks.cmp.pd=sweep{jj}.pks.cmp.pd(min_dist_cmp_to_full>peak_distance_thresh_cmp_full_min);
%         sweep{jj}.pks.all.pzt=[sweep{jj}.pks.cmp.pzt,sweep{jj}.pks.full.pzt];
%         sweep{jj}.pks.all.pd=[sweep{jj}.pks.cmp.pd,sweep{jj}.pks.full.pd];
%         sweep{jj}.pks.all.min_pzt_v_diff=min_diff(sweep{jj}.pks.all.pzt);
% 
% 
% %         if args.plot.all || (args.plot.failed && ~single_mode_vec(jj))
%             title_strs={'Failed','OK'};
%             sfigure(1);
%             subplot(2,2,1)
%             cla;
%              set(gcf,'color','w')
%             title('smoothing pzt v')
%             time=(1:numel(sub_ptz_raw))/sr;
%             plot(time*1e3,sub_ptz_raw,'r')
%             hold on
%             plot(time*1e3,sub_dat_smooth,'k')
%             hold off
%             xlabel('time (ms)')
%             ylabel('volts (v)')
% 
%             subplot(2,2,2)
%             cla;
%             %select a single scan
%             %to do this we will mask out the positive slope
%             diff_sub_raw=diff(sub_ptz_raw)*sr;
%             diff_time=time(1:end-1)+0.5*(time(2)-time(1));
%             plot(diff_time*1e3,diff_sub_raw)
%             hold on
%             plot(diff_time*1e3,sub_dat_grad_smooth)
%             plot(diff_time*1e3,pos_slope_mask*2e3)
%             hold off
%             ylim([min(diff_sub_raw),max(diff_sub_raw)])
%             xlabel('time (ms)')
%             ylabel('grad (v/s)')
% 
% 
%             subplot(2,2,3)
%             cla;
%             plot(sweep{jj}.pzt_smooth,sweep{jj}.pd_full_raw,'k')
%             xlabel('pzt(v)')
%             ylabel('pd (v)')
%             hold on
%             plot(sweep{jj}.pks.full.pzt,sweep{jj}.pks.full.pd,'xk','markersize',20)
% %             plot(sweep{jj}.pzt_smooth,sweep{jj}.pd_cmp_raw*args.cmp_multiplier_disp,'r')
% %             plot(sweep{jj}.pks.cmp.pzt,sweep{jj}.pks.cmp.pd*args.cmp_multiplier_disp,'rx','markersize',20);
%             hold off
%             
% 
%             pause(1e-6)
% %         end %end plots
%     end %enough points to work with
% 
% if jj>=jjmax || ~single_mode_vec(jj) %if something is not single mode do not continue
%     test_sm_while=false;
% end 
% end%end while loop over checking mode
% if sum(~single_mode_vec)==0
%     ai_log_single_out.single_mode=true;
% else
%     ai_log_single_out.single_mode=false;
% end

    
    