%%%%%%%%%%%%%%%%%%%%%%%%%%
% A script to analyse the data we collected from a pulsed atom laser which measures how the trap freq changes with
% application of the tune out probe beam.
% The script:
%   * Imports all the data files into a data structure
%   * Imports the log file, which indexes all the data as a) data run,
%       including wavelength, or b) calibration run to check the control
%       signal
%   * extracts the response at each set point by fitting the aliased trap freq
%   * Produces a plot of the measured response vs probe
%       wavelength, identifying the tuneout wavelength and giving a
%       stastistical uncertainty.
%the data structure
%   first level is instrument or anal method
%   try and pass things between the modules of code only using the 'data' struct


% Other m-files required: import_data,find_data_files,dld_raw_to_txy,masktxy,data_tcreate,
%                         dld_read_5channels_reconst_multi_imp,constants
% Also See:
% Subfunctions: none
% MAT-files required: none
%
% Known BUGS/ Possible Improvements
%   -the fit error depends on wavelength indicating that the model does not
%   have enough freedom
%   -save analysis results
%	-make plots more compact
%   -harmonize the anal opts
%	-place more sections into functions
%	-clean up the fit section
%	-write a n depth function cashing wrapper with hash lookup
%		-alow partial updates
%       
%
% Author: Bryce Henson
% email: Bryce.Henson@live.com
% Last revision:2018-10-01


%close all
%clear all
tic
%%
% BEGIN USER VAR-------------------------------------------------
%tdc_import_opts.dir='Y:\EXPERIMENT-DATA\Tune Out V2\20180826_testing_wm_log\';
%tdc_import_opts.dir='\\amplpc29\Users\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output';
%tdc_import_opts.dir='\\amplpc29\Users\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20180829_half_wp_353';
tdc_import_opts.dir='Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20181002_halfwp_236_stab3\';
tdc_import_opts.file_name='d';
tdc_import_opts.force_load_save=false;   %takes precidence over force_reimport
tdc_import_opts.force_reimport=false;
tdc_import_opts.force_forc=false;
tdc_import_opts.dld_xy_rot=0.61;
%Should probably try optimizing these
tmp_xlim=[-30e-3, 30e-3];     %tight XY lims to eliminate hot spot from destroying pulse widths
tmp_ylim=[-30e-3, 30e-3];
tlim=[0,4];
tdc_import_opts.txylim=[tlim;tmp_xlim;tmp_ylim];

anal_opts.max_runtime=inf;%cut off the data run after some number of hours
anal_opts.atom_laser.pulsedt=8.000e-3;
anal_opts.atom_laser.t0=0.41784; %center i ntime of the first pulse
anal_opts.atom_laser.start_pulse=1; %atom laser pulse to start with
anal_opts.atom_laser.pulses=200;
anal_opts.atom_laser.appr_osc_freq_guess=[52,40,40];
anal_opts.atom_laser.pulse_twindow=anal_opts.atom_laser.pulsedt*0.9;
anal_opts.atom_laser.xylim=tdc_import_opts.txylim(2:3,:); %set same lims for pulses as import
anal_opts.fall_time=0.417;
anal_opts.qe=0.09;

anal_opts.trig_dld=20.3;
anal_opts.dld_aquire=4;
anal_opts.trig_ai_in=20;


anal_opts.wm_log.plot_all=true;
anal_opts.wm_log.plot_failed=true;

histplot.binsx=1000;
histplot.blur=1;
histplot.xlim=[-20,20]*1e-3;
histplot.tlim=[0.86,1.08];
histplot.dimesion=2; %Select coordinate to bin. 1=X, 2=Y.

% END USER VAR-----------------------------------------------------------
data=[]; %CLEAR THE DATA
anal_out=[];
%set upa an output dir %https://gist.github.com/ferryzhou/2269380
if tdc_import_opts.dir(end) ~= '\', dirpath = [dirpath '\']; end
if (exist([tdc_import_opts.dir,'out'], 'dir') == 0), mkdir([tdc_import_opts.dir,'out']); end
 
anal_out.dir=sprintf('%sout\\%s\\',...
    tdc_import_opts.dir,datestr(datetime('now'),'yyyymmddTHHMMSS'));
if (exist(anal_out.dir, 'dir') == 0), mkdir(anal_out.dir); end
 
diary([anal_out.dir,'anal.txt'])

%sets up the struct 'data' which will contain everything you could want incuding the txy data and
%the information from the logs
addpath('Colormaps') 
addpath('FileTime_29Jun2011') %used for high precision windows timestamps in import_data
constants
anal_opts.velocity=const.g0*anal_opts.fall_time;
tdc_import_opts.shot_num=find_data_files(tdc_import_opts);
[mcp_tdc_data,import_opts]=import_mcp_tdc_data(tdc_import_opts);
data.mcp_tdc=mcp_tdc_data;

%TO DO FUNCTIONALIZE
%import the wavemeter log
%adaptively to deal with the 2 different log files that are in the data
clear('wm_log')
wm_log.dir = strcat(tdc_import_opts.dir,'log_LabviewMatlab.txt');
fid = fopen(wm_log.dir );
wm_log.cell=textscan(fid,'%s','Delimiter','\n');
fclose(fid);
wm_log.cell=wm_log.cell{1};
for ii=1:size(wm_log.cell,1)
    if ~isequal(wm_log.cell{ii},'') %catch the empty case
        if contains(wm_log.cell{ii},'measure_probe')
            line_cells=textscan(wm_log.cell{ii},'%f %s %s %s %f %s %u','Delimiter',',');
            wm_log.setpoints(ii)=line_cells{5};
            wm_log.probe_calibration(ii)=false;
            wm_log.iter_nums(ii)=line_cells{7};
        elseif contains(wm_log.cell{ii},'calibrate')
            line_cells=textscan(wm_log.cell{ii},'%f %s %s %s %s %u','Delimiter',',');
            wm_log.setpoints(ii)=NaN;
            wm_log.probe_calibration(ii)=true;
            wm_log.iter_nums(ii)=line_cells{6};
        else %deals with the legacy case (only 20180813_CW_AL_tuneout_scan)
            line_cells=textscan(wm_log.cell{ii},'%f %s %s %s %f %s %u','Delimiter',',');
            wm_log.setpoints(ii)=line_cells{5};
            wm_log.probe_calibration(ii)=false;
            wm_log.iter_nums(ii)=line_cells{7};
        end
        wm_log.posix_times(ii)=line_cells{1};
        wm_log.iso_times{ii}=line_cells{2};
    end
end
data.labview.setpoint=wm_log.setpoints*1e6; %convert to hz
data.labview.time=wm_log.posix_times;
data.labview.shot_num=wm_log.iter_nums;
data.labview.calibration=wm_log.probe_calibration;
%can check that the times look ok
% plot(data.mcp_tdc.write_time-data.probe.time)

%%
%try to import the wm log files
%should ckeck that the same files are avail when import from save
wm_log_import_opts.dir=tdc_import_opts.dir;
wm_log_import_opts.force_reimport=false;
wm_log_name='log_wm_';
wm_logs=dir([wm_log_import_opts.dir,wm_log_name,'*.txt']);
wm_log_import_opts.names={wm_logs.name};
data.wm_log=wm_log_import(wm_log_import_opts);

%%
sfigure(1)
%create a list of indicies (of the mcp_tdc) that have an ok number of counts
%exclude the very low and then set the thresh based on the sd of the remaining
not_zero_files=data.mcp_tdc.num_counts>1e3; 
num_thresh=mean(data.mcp_tdc.num_counts(not_zero_files))-4*std(data.mcp_tdc.num_counts(not_zero_files));
data.mcp_tdc.num_ok=data.mcp_tdc.num_counts>num_thresh & ...
    (data.mcp_tdc.time_create_write(:,1)'-data.mcp_tdc.time_create_write(1,1))<(anal_opts.max_runtime*60*60);
fprintf('shots number ok %u out of %u \n',sum(data.mcp_tdc.num_ok),numel(data.mcp_tdc.num_ok))

plot((data.mcp_tdc.time_create_write(:,2)-data.mcp_tdc.time_create_write(1,2))/(60*60),data.mcp_tdc.num_counts)
xlabel('time (h)')
ylabel('total counts')
title('num count run trend')

%% Match up shot numbers
%because the mcp-dld detector has no direct communication with the bec computer
% the data.labview.shot_num does not nessesarily correspond to data.mcp_tdc.shot_num

match_times=true;
if match_times
    %try and match up the file with if it is a calibaration using the time
    %it is slightly overkill here to search each one, but just being extra
    %cautious/flexible
    time_thresh=2; %how close for the times to be considered the same shot
    %lets examine what the time difference does
    sfigure(45);
    set(gcf,'color','w')
    clf
    imax=min([size(data.labview.time,2),size(data.mcp_tdc.time_create_write,1)]);
    %imax=5000;
    plot(data.mcp_tdc.time_create_write(1:imax,1)'-data.labview.time(1:imax)...
        -anal_opts.trig_dld)
    mean_delay_labview_tdc=mean(data.labview.time(1:imax)-data.mcp_tdc.time_create_write(1:imax,1)');
    xlabel('shot number')
    ylabel('time between labview and mcp tdc')
    %to do include ai_log
    iimax=size(data.mcp_tdc.time_create_write(:,1),1);
    data.mcp_tdc.probe.calibration=nan(iimax,1);
    data.mcp_tdc.labview_shot_num=nan(iimax,1);
    %loop over all the tdc_files 
    for ii=1:iimax
        %predict the labview master trig time
        %use the write time to handle being unpacked from 7z
        est_labview_start=data.mcp_tdc.time_create_write(ii,2)...
            -anal_opts.trig_dld-anal_opts.dld_aquire;
        [tval,nearest_idx]=closest_value(data.labview.time...
            ,est_labview_start);
        if abs(tval-est_labview_start)<time_thresh
            data.mcp_tdc.labview_shot_num(ii)=data.labview.shot_num(nearest_idx);
            data.mcp_tdc.probe.calibration(ii)=data.labview.calibration(nearest_idx);
        end 
    end
else
    %just do it the boring way and hope that the tdc was set up right and
    %there were no false trigers
    %TO DO
    
end

%% use the analog input log to check if the probe beam pd voltage was ok and
%that the laser was single mode
anal_opts.ai_log.dir=tdc_import_opts.dir;
anal_opts.ai_log.force_reimport=false ;
anal_opts.ai_log.force_load_save=false;
anal_opts.ai_log.log_name='log_analog_in_';
anal_opts.ai_log.pd.set=2;
anal_opts.ai_log.pd.diff_thresh=2;
anal_opts.ai_log.pd.std_thresh=0.1;
anal_opts.ai_log.pd.time_start=0.2;
anal_opts.ai_log.pd.time_stop=2;
anal_opts.ai_log.sfp.num_checks=10; %how many places to check that the laser is single mode
anal_opts.ai_log.sfp.thresh_cmp_peak=20e-3; %theshold on the compressed signal to be considered a peak
anal_opts.ai_log.sfp.peak_dist_min_pass=4.5;%minimum (min difference)between peaks for the laser to be considered single mode
anal_opts.ai_log.plot.all=false;
anal_opts.ai_log.plot.failed=true;


%because im only passing the ai_log feild to aviod conflicts forcing a reimport i need to coppy these feilds
anal_opts.ai_log.trig_dld=anal_opts.trig_dld;
anal_opts.ai_log.trig_ai_in=anal_opts.trig_ai_in;
[data,ai_log_out]=ai_log_import(anal_opts.ai_log,data);
%copy the output across
data.ai_log=ai_log_out.ai_log;
data.mcp_tdc.probe.ok.reg_pd=ai_log_out.mcp_tdc.probe.ok.reg_pd;
data.mcp_tdc.probe.ok.sfp=ai_log_out.mcp_tdc.probe.ok.sfp;


%save([datestr(datetime('now'),'yyyymmddTHHMMSS'),'.mat'],'-v7.3')


%%
%so we want to define a time window
%-for checking if the doubler was ok
%averaging the wavelength of the laser
%i think it will be anal_opts.atom_laser.t0 after the creation time of the tdc file
%compexity is that the time that the tdc file is wrote/reated is not relaible and depend on the flux rate and avaialble mem
%to this end find the closest labview update time and go back one then fowards

anal_opts.wm_log.plot_all=false;
anal_opts.wm_log.plot_failed=true;

time_pd_padding=4; %check this many s each side of probe
time_blue_padding=1; %check this many seconde each side of probe
time_probe=3;
ecd_volt_thresh=0.5;

red_sd_thresh=12;
red_range_thresh=50;
rvb_thresh=10;

mean_shot_duration=mean(diff(data.labview.time(:)));
data.mcp_tdc.probe_freq=[];
fprintf('finding mean wavelengths for files %04u:%04u',size(data.mcp_tdc.time_create_write(:,1),1),0)
if ~issorted(data.wm_log.feedback.posix_time) || ~issorted(data.wm_log.read_all_adc.posix_time)
    error('binary search wont work')
end
iimax=size(data.mcp_tdc.time_create_write(:,1),1);

data.mcp_tdc.probe.freq.set=nan(iimax,1);
data.mcp_tdc.probe.freq.act.mean=nan(iimax,1);
data.mcp_tdc.probe.freq.act.std=nan(iimax,1);
data.mcp_tdc.probe.freq.rvb.mean=nan(iimax,1);

%GUILTY UNTILL PROVEN INOCENT  !!!!!!!!!!!!!!!!!
data.mcp_tdc.probe.ok.freq=false(iimax,1); %frequency reading
data.mcp_tdc.probe.ok.rvb=false(iimax,1);  %2r-b check
data.mcp_tdc.probe.ok.ecd_pd=false(iimax,1);  %ecd pd value
  

sfigure(11);
set(gcf,'color','w')

for ii=1:iimax
    est_labview_start=data.mcp_tdc.time_create_write(ii,2)...
            -anal_opts.trig_dld-anal_opts.dld_aquire;
    %find the labview update that was nearest to the tdc_time minus the tdc_trig time
    %this should be the labview itteration that made this data 
    %is this the best way to do it?
    [time_nearest_lv,idx_nearest_lv]=closest_value(data.labview.time,est_labview_start);
    abs_anal_opts.trig_dld_on_main_comp=data.labview.time(idx_nearest_lv)+anal_opts.trig_dld;
    time_lower=abs_anal_opts.trig_dld_on_main_comp+anal_opts.atom_laser.t0-anal_opts.fall_time-0.5; %the probe turns on
    time_upper=abs_anal_opts.trig_dld_on_main_comp+anal_opts.atom_laser.t0-anal_opts.fall_time+time_probe+0.5; %when the probe beam goes off
    
    %CHECK PD VALUE
    %check if there were some readings of the ecd output voltage withing time_padding seconds of this period
    %the padding is because EDC readings are only taken every 3s
    ecd_pd_mask_idx=fast_sorted_mask(data.wm_log.read_all_adc.posix_time,...
                time_lower-time_pd_padding,...
                time_upper+time_pd_padding);
    ecd_pd_measurments=max(0,ecd_pd_mask_idx(2)-ecd_pd_mask_idx(1));
    pd_readings=data.wm_log.read_all_adc.value10(ecd_pd_mask_idx(1):ecd_pd_mask_idx(2));  
    if sum(ecd_pd_measurments)>=2 && sum(pd_readings<ecd_volt_thresh)==0 %if there are not a few during this time then throw an error
        data.mcp_tdc.probe.ok.ecd_pd(ii)=true;
    end
    
    %CHECK wm freq
    %could check if close to set pt, but for now will just check if flat
    wm_red_mask_idx=fast_sorted_mask(data.wm_log.feedback.posix_time,time_lower,time_upper);
    red_freqs=data.wm_log.feedback.actual(wm_red_mask_idx(1):wm_red_mask_idx(2));
    data.mcp_tdc.probe.freq.set(ii)=mean(data.wm_log.feedback.setpt(wm_red_mask_idx(1):wm_red_mask_idx(2)));
    data.mcp_tdc.probe.freq.act.mean(ii)=mean(red_freqs);
    data.mcp_tdc.probe.freq.act.std(ii)=std(red_freqs);
    data.mcp_tdc.probe.freq.error(ii)=false;
    %             break
    if numel(red_freqs)~=0 && std(red_freqs)<red_sd_thresh && range(red_freqs)<red_range_thresh
       data.mcp_tdc.probe.ok.freq(ii)=true;
    end
    
    %check that the blue is ~ 2x the red freq
    wm_blue_mask_idx=fast_sorted_mask(data.wm_log.blue_freq.posix_time,...
                time_lower-time_blue_padding,...
                time_upper+time_blue_padding);
    blue_freqs=data.wm_log.blue_freq.value(wm_blue_mask_idx(1):wm_blue_mask_idx(2));
    
    data.mcp_tdc.probe.freq.rvb.mean(ii)=mean(blue_freqs)-data.mcp_tdc.probe.freq.act.mean(ii)*2;
    if abs(data.mcp_tdc.probe.freq.rvb.mean(ii))<rvb_thresh
        data.mcp_tdc.probe.ok.rvb(ii)=true;
    end
    %do some plots if plot_all or if plot_failed what failed
    if anal_opts.wm_log.plot_all || (anal_opts.wm_log.plot_failed &&...
       (~data.mcp_tdc.probe.ok.freq(ii) || ~data.mcp_tdc.probe.ok.rvb(ii) ||...
       ~data.mcp_tdc.probe.ok.ecd_pd(ii)))
        
        sfigure(11);
        subplot(3,1,1)
        plot(data.wm_log.read_all_adc.posix_time(ecd_pd_mask_idx(1):...
            ecd_pd_mask_idx(2))+anal_opts.atom_laser.t0-...
            abs_anal_opts.trig_dld_on_main_comp,pd_readings)
        xlabel('Probe on time(s)')
        ylabel('Doubler PD voltage')
        title('Doubler output check')
        subplot(3,1,2)
        tmp_wav_avg=mean(red_freqs);
        tmp_wav_cen=red_freqs-tmp_wav_avg;
        plot(anal_opts.atom_laser.t0+ data.wm_log.feedback.posix_time(wm_red_mask_idx(1):wm_red_mask_idx(2))...
            -abs_anal_opts.trig_dld_on_main_comp,tmp_wav_cen)
        title(sprintf('Red freq - %.1fMHz',tmp_wav_avg))
        xlabel('Probe on time(s)')
        ylabel('Red Freq (MHz)')
        yl=ylim;

        
        subplot(3,1,3)
        %should really interpolate here to give the difference properly as
        %a function of time
        title('2red(mean)-blue difference')
        tmp_rvb=blue_freqs-data.mcp_tdc.probe.freq.act.mean(ii)*2;
        plot(anal_opts.atom_laser.t0+ data.wm_log.blue_freq.posix_time(wm_blue_mask_idx(1):wm_blue_mask_idx(2))...
            -abs_anal_opts.trig_dld_on_main_comp,tmp_rvb)
        xlabel('Probe on time(s)')
        ylabel('2r-b Freq (MHz)')
        title('Blue red difference')
        pause(1e-2)
    end
    if mod(ii,1e2)==0,fprintf('\b\b\b\b%04u',ii),end
end
fprintf('...Done\n')


%%
%Here we will do a plot of all the checks and then combine them into one
%master 'ok'/check vector
%data.mcp_tdc.probe.ok.reg_pd
%data.mcp_tdc.probe.ok.sfp
%data.mcp_tdc.num_ok

%data.mcp_tdc.probe.ok.freq; %frequency reading is stable
%data.mcp_tdc.probe.ok.rvb;  %2r-b check
%data.mcp_tdc.probe.ok.ecd_pd;  %ecd pd value
  

figure(12);
clf
set(gcf,'color','w')
subplot(2,1,1)
%plot all the logics, dither it a bit to make it easier to figure out
%culprits
line_width=1.5;
stairs(data.mcp_tdc.shot_num,data.mcp_tdc.num_ok-0.03,'LineWidth',line_width)
hold on
stairs(data.mcp_tdc.shot_num,data.mcp_tdc.probe.ok.reg_pd-0.01,'LineWidth',line_width)
stairs(data.mcp_tdc.shot_num,data.mcp_tdc.probe.ok.sfp-0.02,'LineWidth',line_width)
stairs(data.mcp_tdc.shot_num,data.mcp_tdc.probe.ok.freq+0.00,'LineWidth',line_width)
%forms another check that the laser is single mode
stairs(data.mcp_tdc.shot_num,data.mcp_tdc.probe.ok.rvb+0.01,'LineWidth',line_width)
 %this is reduncant as it should be caught by the reg_pd measurement
stairs(data.mcp_tdc.shot_num,data.mcp_tdc.probe.ok.ecd_pd+0.02,'LineWidth',line_width) 
tmp_cal=data.mcp_tdc.probe.calibration;
tmp_cal(~isnan(tmp_cal))=~tmp_cal(~isnan(tmp_cal));
stairs(data.mcp_tdc.shot_num,tmp_cal-0.02,'LineWidth',line_width) 
hold off
title('Checks')
xlabel('Shot Number')
ylabel('Good?')
set(gca,'ytick',[0,1],'yticklabel',{'False','True'})
ylim([-0.1,1.1])
xl=xlim;
xlim([1,xl(2)])
legend('number','pd reg','single mode & pd','freq stable','RvB','ecd pd(ignored)','NOT(calibration)')
yticks([0 1])
subplot(2,1,2)
tmp_cal=data.mcp_tdc.probe.calibration;
tmp_cal(isnan(tmp_cal))=false;
%must have good atom number AND (good probe OR be calibration)
tmp_probe_ok=(data.mcp_tdc.probe.ok.sfp &...
    data.mcp_tdc.probe.ok.reg_pd &...
    data.mcp_tdc.probe.ok.freq &....
    data.mcp_tdc.probe.ok.rvb );
tmp_all_ok=data.mcp_tdc.num_ok' &...
    (tmp_probe_ok| tmp_cal);
    %data.mcp_tdc.probe.ok.ecd_pd;
stairs(data.mcp_tdc.shot_num,tmp_all_ok,'LineWidth',line_width)
hold on
stairs(data.mcp_tdc.shot_num,tmp_probe_ok+0.01,'LineWidth',line_width)
hold off
legend('all','probe')
ylabel('Good?')
xlabel('Shot Number')
set(gca,'ytick',[0,1],'yticklabel',{'False','True'})
title('ALL ok')
ylim([-0.1,1.1])
tmp_num_shots=numel(data.mcp_tdc.shot_num);
tmp_num_ok_shots=sum(tmp_all_ok);
fprintf('ok logic gives %u / %u shots for yeild %04.1f %%\n',...
    tmp_num_ok_shots,tmp_num_shots,1e2*tmp_num_ok_shots/tmp_num_shots)

set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1600, 900])
data.mcp_tdc.all_ok=tmp_all_ok;
plot_name='check_logics';
saveas(gcf,[anal_out.dir,plot_name,'.png'])
saveas(gcf,[anal_out.dir,plot_name,'.fig'])

%% binning the pulses
%first step is to bin up each pulse of the atom laser

tic
iimax=size(data.mcp_tdc.counts_txy,2);
data.mcp_tdc.al_pulses=[];
data.mcp_tdc.al_pulses.pulsedt=anal_opts.atom_laser.pulsedt;
data.mcp_tdc.al_pulses.window=nan(anal_opts.atom_laser.pulses,3,2); %initalize
data.mcp_tdc.al_pulses.num_counts=nan(iimax,anal_opts.atom_laser.pulses);
  
plots=false;
fprintf('binning pulses in files %04u:%04u',size(data.mcp_tdc.counts_txy,2),0)
first_good_shot=true;
for shot=1:iimax
        if data.mcp_tdc.all_ok(shot)
            for pulse=1:anal_opts.atom_laser.pulses
                %set up time window centered arround t0
                trange=anal_opts.atom_laser.t0+anal_opts.atom_laser.pulsedt...
                    *(anal_opts.atom_laser.start_pulse+pulse-2)+...
                    anal_opts.atom_laser.pulse_twindow*[-0.5,0.5];
                pulse_win_txy=[trange;anal_opts.atom_laser.xylim]; 
                counts_pulse=masktxy(data.mcp_tdc.counts_txy{shot},pulse_win_txy);
                if plots
                    sfigure(79);
                    set(gcf,'Color',[1 1 1]);
                    subplot(3,1,1)
                    hist(counts_pulse(:,1),100)
                    xlabel('t')
                    title('full')
                    subplot(3,1,2)
                    hist(counts_pulse(:,2),100)
                    xlabel('x')
                    title('full')
                    subplot(3,1,3)
                    hist(counts_pulse(:,3),100)
                    xlabel('y')
                    title('full')
                    pause(0.01)
                end
                if first_good_shot
                    
                    %only need to store this on first shot becasue the same for
                    %all shots
                    data.mcp_tdc.al_pulses.window(pulse,:,:)=pulse_win_txy; 
                    data.mcp_tdc.al_pulses.time(pulse,:)=(pulse+anal_opts.atom_laser.start_pulse-2)*anal_opts.atom_laser.pulsedt;
                end
                data.mcp_tdc.al_pulses.num_counts(shot,pulse)=size(counts_pulse(:,3),1);
                data.mcp_tdc.al_pulses.pos_stat(shot,pulse,:)=[...
                                           mean(counts_pulse(:,1)),...
                                           mean(counts_pulse(:,2)),...
                                           mean(counts_pulse(:,3)),...
                                           std(counts_pulse(:,1)),...
                                           std(counts_pulse(:,2)),...
                                           std(counts_pulse(:,3))]; 
            end%pulse
        if first_good_shot,first_good_shot=false; end
        end%is data.mcp_tdc.all_ok
        if mod(shot,10)==0,fprintf('\b\b\b\b%04u',shot),end    
%to set the pulse t0 right it can be handy to uncomment the next line
%fprintf('\nmean time %3.5f            \n ',mean(data.mcp_tdc.al_pulses.pos_stat(shot,:,1)-data.mcp_tdc.al_pulses.time(:)'))
end%shots
fprintf('...Done\n') 

toc



%% Fitting the Trap freq
%using the binned data we fit the trap freq to each shot
%loop over every shot in data.mcp_tdc but output nothing if
%data.mcp_tdc.all_ok(ii)=false
%for compactness could use max((1:numel(data.mcp_tdc.all_ok)).*data.mcp_tdc.all_ok');
%find the last good shot, but this fuck up the mask && mask operations
%later
iimax=size(data.mcp_tdc.counts_txy,2); 
plots=false;
anal_opts.atom_laser.appr_osc_freq_guess=[52,46.7,40];
 %try and guess the trap freq based on the peak of the fft, needed when the
 %kick amp is changing
adaptive_fit_freq=true;
%ignore some of the fit errors
warning('off','stats:nlinfit:ModelConstantWRTParam');
warning('off','MATLAB:rankDeficientMatrix');
data.osc_fit=[]; %clear the output struct
%prealocate so that i can do the logics later
data.osc_fit.dld_shot_idx=nan(1,iimax);
data.osc_fit.model_coefs=nan(iimax,8,2);
data.osc_fit.fit_rmse=nan(1,iimax);
data.osc_fit.model=cell(1,iimax);
fprintf('Fitting oscillations in shots %04i:%04i',iimax,0)
for ii=1:iimax
    %position that data appears in data.mcp_tdc, not ness shot number
    %specified because we will remove any elements of osc_fit that did not
    %fit because of all_ok condition
    data.osc_fit.dld_shot_idx(ii)=ii;
    %shot number eg d123.txt as recorded by the tdc computer, not ness lv
    %number
    dld_shot_num=data.mcp_tdc.shot_num(ii);
    if data.mcp_tdc.all_ok(ii)
        data.osc_fit.dld_shot_num(ii)=dld_shot_num;
        %construct a more convinent temp variable txyz_tmp wich is the position in mm for use in the fit
        x_tmp=1e3*squeeze(data.mcp_tdc.al_pulses.pos_stat(ii,:,2));
        x_tmp=x_tmp-nanmean(x_tmp);
        y_tmp=1e3*squeeze(data.mcp_tdc.al_pulses.pos_stat(ii,:,3));
        y_tmp=y_tmp-nanmean(y_tmp);
        z_tmp=data.mcp_tdc.al_pulses.time-squeeze(data.mcp_tdc.al_pulses.pos_stat(ii,:,1))';
        z_tmp=z_tmp'*anal_opts.velocity*1e3;
        z_tmp=z_tmp-nanmean(z_tmp);
        txyz_tmp=[data.mcp_tdc.al_pulses.time';x_tmp;y_tmp;z_tmp];
        sqrtn=sqrt(data.mcp_tdc.al_pulses.num_counts(ii,:)); %find the statistical uncert in a single shot
        xerr=1e3*squeeze(data.mcp_tdc.al_pulses.pos_stat(ii,:,5))./sqrtn;
        yerr=1e3*squeeze(data.mcp_tdc.al_pulses.pos_stat(ii,:,6))./sqrtn;
        zerr=1e3*squeeze(data.mcp_tdc.al_pulses.pos_stat(ii,:,4))*anal_opts.velocity./sqrtn;
        xyzerr_tmp=[xerr;yerr;zerr];
        xyzerr_tmp(:,sqrtn<2)=nan;

        %remove any data pts with nan position
        mask=sum(isnan(txyz_tmp),1)==0;
        xyzerr_tmp=xyzerr_tmp(:,mask);
        txyz_tmp=txyz_tmp(:,mask);

        %try to find the peak osc freq to start the fit there
        if adaptive_fit_freq
            out=fft_tx(txyz_tmp(1,:),txyz_tmp(histplot.dimesion+1,:),10);
            [~,nearest_idx]=max(abs(out(2,:)));
            fit_freq=out(1,nearest_idx);
            %fft_phase=angle(out(2,nearest_idx))+0.535;
        else
            fit_freq=anal_opts.atom_laser.appr_osc_freq_guess(histplot.dimesion);
        end
        modelfun = @(b,x) exp(-x(:,1).*max(0,b(7))).*b(1).*sin(b(2)*x(:,1)*pi*2+b(3)*pi*2)+b(4)+b(8)*x(:,1)+b(5)*x(:,2)+b(6)*x(:,3);
        beta0=[std(txyz_tmp(histplot.dimesion+1,:))*8, fit_freq, 0, 1,0,0,2,0.01];
        cof_names={'amp','freq','phase','offset','ycpl','zcpl','damp','grad'};
        opt = statset('TolFun',1e-10,'TolX',1e-10,'MaxIter',1e4,...
            'UseParallel',1);
        %select the aproapriate values to go in the response variable
        idx=1:4;
        idx(histplot.dimesion+1)=[];
        predictor=txyz_tmp(idx,:)';
        weights=1./(xyzerr_tmp(histplot.dimesion,:).^2);
        weights(isnan(weights))=1e-20; %if nan then set to smallest value you can
        weights=weights/sum(weights);
        %predictor=[tvalues,xvalues,zvalues];
        fitobject=fitnlm(predictor,txyz_tmp(histplot.dimesion+1,:)',modelfun,beta0,...
            'Weights',weights,'options',opt,...
            'CoefficientNames',cof_names);
        data.osc_fit.model{ii}=fitobject;
        fitparam=fitobject.Coefficients;
        data.osc_fit.model_coefs(ii,:,:)=[fitparam.Estimate,fitparam.SE];
        data.osc_fit.fit_rmse(ii)=fitobject.RMSE;
        %limiting frequnecy prediction from http://adsabs.harvard.edu/full/1999DSSN...13...28M
        meanwidth=sqrt(mean(squeeze(data.mcp_tdc.al_pulses.pos_stat(ii,:,5)).^2))*1e3;
        frequnclim=sqrt(6/sum(data.mcp_tdc.al_pulses.num_counts(ii,:)))*...
            (1/(pi*range(data.mcp_tdc.al_pulses.time)))*...
            (meanwidth/fitparam{2,1});
        %fprintf('sampling limit %2.3g Hz, fint unc %2.3g Hz, ratio %2.3g \n',[frequnclim,fitparam{2,2},fitparam{2,2}/frequnclim])
        data.osc_fit.fit_sample_limit{ii}=[frequnclim,fitparam{2,2},fitparam{2,2}/frequnclim];
        if plots
            tplotvalues=linspace(min(data.mcp_tdc.al_pulses.time),...
                max(data.mcp_tdc.al_pulses.time),1e5)';
            predictorplot=[tplotvalues,...
                       interp1(predictor(:,1),predictor(:,2),tplotvalues),...
                       interp1(predictor(:,1),predictor(:,3),tplotvalues)];
            [prediction,ci]=predict(fitobject,predictorplot);
            sfigure(2);
            clf
            set(gcf,'color','w')
            subplot(2,1,1)
            plot(txyz_tmp(1,:),txyz_tmp(2,:),'kx-')
            hold on
            plot(txyz_tmp(1,:),txyz_tmp(3,:),'rx-')
            plot(txyz_tmp(1,:),txyz_tmp(4,:),'bx-')
            hold off
            ylabel('X Pos (mm)')
            xlabel('Time (s)')
            set(gca,'Ydir','normal')
            set(gcf,'Color',[1 1 1]);
            legend('x','y','z')
            pause(0.05)

            subplot(2,1,2)
            plot(predictorplot(:,1),prediction,'-','LineWidth',1.5,'Color',[0.5 0.5 0.5])
            ax = gca;
            set(ax, {'XColor', 'YColor'}, {'k', 'k'});
            hold on
            plot(predictorplot(:,1),ci(:,1),'-','LineWidth',1.5,'Color','k')
            plot(predictorplot(:,1),ci(:,2),'-','LineWidth',1.5,'Color','k')
            errorbar(predictor(:,1),txyz_tmp(histplot.dimesion+1,:)',xyzerr_tmp(histplot.dimesion,:),'k.','MarkerSize',10,'CapSize',0,'LineWidth',1,'Color','r') 
            set(gcf,'Color',[1 1 1]);
            ylabel('X(mm)')
            xlabel('Time (s)')
            hold off
            ax = gca;
            set(ax, {'XColor', 'YColor'}, {'k', 'k'});
            set(gca,'linewidth',1.0)
            saveas(gca,sprintf('.\\out\\fit_dld_shot_num%04u.png',dld_shot_num))
        end% PLOTS
    end
    fprintf('\b\b\b\b%04u',ii)
end
fprintf('...Done\n')

data.osc_fit.ok.did_fits=~cellfun(@(x) isequal(x,[]),data.osc_fit.model);
        

%%
%look for failed fits
fprintf('mean fit error %f\n',...
    mean(data.osc_fit.fit_rmse(data.osc_fit.ok.did_fits)))
figure(1)
clf
subplot(2,1,1)
hist(data.osc_fit.fit_rmse(data.osc_fit.ok.did_fits))
xlabel('RMSE')
ylabel('counts')
subplot(2,1,2)
plot(data.osc_fit.fit_rmse(data.osc_fit.ok.did_fits))
xlabel('shot idx')
mask=data.osc_fit.ok.did_fits;
data.osc_fit.ok.rmse=mask & data.osc_fit.fit_rmse...
    < nanmean(data.osc_fit.fit_rmse(mask))+2*nanstd(data.osc_fit.fit_rmse(mask));
data.osc_fit.ok.rmse(data.osc_fit.ok.rmse)=abs(data.osc_fit.model_coefs(data.osc_fit.ok.rmse,2,1)'...
    -mean(data.osc_fit.model_coefs(data.osc_fit.ok.rmse,2,1)))<1;
%% Investigate Spurrious fit correlations
sfigure(852);
clf
set(gcf,'color','w')
%see if the fit error depends on the probe freq
subplot(3,3,1)
tmp_probe_freq=data.mcp_tdc.probe.freq.act.mean(data.osc_fit.ok.rmse);
tmp_probe_freq=(tmp_probe_freq-nanmean(tmp_probe_freq))*1e-3;
plot(tmp_probe_freq,...
    data.osc_fit.fit_rmse(data.osc_fit.ok.rmse),'xk')
xlabel('probe beam freq (GHz)')
ylabel('fit error')
title('Fit Error')
%see if the damping changes
subplot(3,3,2)
plot(tmp_probe_freq,...
   data.osc_fit.model_coefs((data.osc_fit.ok.rmse),1,1),'xk')
xlabel('probe beam freq (GHz)')
ylabel('Osc amp (mm)')
title('Amp')
subplot(3,3,3)
plot(tmp_probe_freq,...
   data.osc_fit.model_coefs((data.osc_fit.ok.rmse),3,1),'xk')
xlabel('probe beam freq (GHz)')
ylabel('Phase (rad)')
title('Phase')
subplot(3,3,4)
plot(tmp_probe_freq,...
   data.osc_fit.model_coefs((data.osc_fit.ok.rmse),4,1),'xk')
xlabel('probe beam freq (GHz)')
ylabel('offset')
title('offset')
subplot(3,3,5)
plot(tmp_probe_freq,...
   data.osc_fit.model_coefs((data.osc_fit.ok.rmse),5,1),'xk')
xlabel('probe beam freq (GHz)')
ylabel('ycpl')
title('ycpl')
subplot(3,3,6)
plot(tmp_probe_freq,...
   data.osc_fit.model_coefs((data.osc_fit.ok.rmse),6,1),'xk')
xlabel('probe beam freq (GHz)')
ylabel('zcpl')
title('zcpl')
%see if the amp changes
subplot(3,3,7)
plot(tmp_probe_freq,...
   1./data.osc_fit.model_coefs((data.osc_fit.ok.rmse),7,1),'xk')
xlabel('probe beam freq (GHz)')
ylabel('Damping Time (s)')
title('Damping')
subplot(3,3,8)
plot(tmp_probe_freq,...
   data.osc_fit.model_coefs((data.osc_fit.ok.rmse),8,1),'xk')
xlabel('probe beam freq (GHz)')
ylabel('Grad mm/s')
title('Grad')
%look for anharminicity with a osc amp/freq correlation
subplot(3,3,9)
plot(abs(data.osc_fit.model_coefs((data.osc_fit.ok.rmse),1,1)),...
    data.osc_fit.model_coefs((data.osc_fit.ok.rmse),2,1),'xk')
xlabel('osc (mm)')
ylabel('fit freq (Hz)')
title('anharmonicity')

set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1600, 900])
data.mcp_tdc.all_ok=tmp_all_ok;
plot_name='fit_correlations';
saveas(gcf,[anal_out.dir,plot_name,'.png'])
saveas(gcf,[anal_out.dir,plot_name,'.fig'])


%% just plot with index
plot(data.osc_fit.model_coefs(data.osc_fit.ok.rmse,2,1)'...
    -mean(data.osc_fit.model_coefs(data.osc_fit.ok.rmse,2,1)))


%% CHECK IF ATOM NUMBER DEPENDS ON PROBE BEAM
figure
clf
plot(data.mcp_tdc.probe.freq.act.mean(data.osc_fit.ok.rmse),data.mcp_tdc.num_counts(data.osc_fit.ok.rmse),'x')


%% plot not calibrations
clf
temp_cal=data.mcp_tdc.probe.calibration';
temp_cal(isnan(temp_cal))=1;
mask=data.osc_fit.ok.rmse & ~temp_cal & ~isnan(data.mcp_tdc.probe.freq.act.mean');
plot(data.mcp_tdc.probe.freq.act.mean(mask),data.osc_fit.model_coefs(mask,2,1),'x')


%% create a model of the underlying trap frequency from the calibrations
addpath('nanconv')
cal_smooth_time=200;
figure(3)
clf
set(gcf,'color','w')
temp_cal=data.mcp_tdc.probe.calibration';
temp_cal(isnan(temp_cal))=0;
cal_dat_mask=data.osc_fit.ok.rmse & temp_cal & ~isnan(data.mcp_tdc.probe.freq.act.mean');
x_tmp=data.mcp_tdc.time_create_write(cal_dat_mask,1);
time_start_cal=x_tmp(1);
x_tmp=x_tmp-time_start_cal;
y_tmp=3*(1/anal_opts.atom_laser.pulsedt)+data.osc_fit.model_coefs(cal_dat_mask,2,1);
freq_drift_model=fit(x_tmp, y_tmp,'poly5');
x_samp=linspace(min(x_tmp),max(x_tmp),1e3);
hour_in_s=1;%60*60;
plot(x_tmp/hour_in_s,y_tmp)
hold on
plot(x_samp/hour_in_s,freq_drift_model(x_samp))
title('Trap Freq Calibration')
xlabel('experiment time (h)')
ylabel('no probe trap freq')
xinterp=linspace(min(x_tmp),max(x_tmp),1e5);
yinterp_raw=interp1(x_tmp,y_tmp,xinterp,'linear');

dx_interp=xinterp(2)-xinterp(1);
kernel = gausswin(ceil(3*cal_smooth_time/(dx_interp)),3);
kernel=kernel/sum(kernel);
yinterp_smooth = nanconv(yinterp_raw,kernel,'edge','1d')';
freq_drift_model=@(x) interp1(xinterp,yinterp_smooth,x-time_start_cal,'linear');
plot(xinterp/hour_in_s,yinterp_raw,'g')
plot(xinterp/hour_in_s,yinterp_smooth,'m')
plot(x_samp/hour_in_s,freq_drift_model(x_samp+time_start_cal),'k')
legend('data','polynomial model','interp data model','smooth data model','freq\_drift\_model' )
hold off

set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1600, 900])
data.mcp_tdc.all_ok=tmp_all_ok;
plot_name='calibration_model';
saveas(gcf,[anal_out.dir,plot_name,'.png'])
saveas(gcf,[anal_out.dir,plot_name,'.fig'])


%%
figure(4)
clf
set(gcf,'color','w')
subplot(2,1,1)
temp_cal=data.mcp_tdc.probe.calibration';
temp_cal(isnan(temp_cal))=1;
probe_dat_mask=data.osc_fit.ok.rmse & ~temp_cal & ~isnan(data.mcp_tdc.probe.freq.act.mean');

probe_freq= data.mcp_tdc.probe.freq.act.mean(probe_dat_mask)*1e6;
corrected_delta=3*(1/anal_opts.atom_laser.pulsedt)+data.osc_fit.model_coefs(probe_dat_mask,2,1)-...
    freq_drift_model(data.mcp_tdc.time_create_write(probe_dat_mask,1));
cdat=viridis(1000);
c_cord=linspace(0,1,size(cdat,1));
shot_time=data.mcp_tdc.time_create_write(probe_dat_mask,1);
shot_time=shot_time-min(shot_time);
shot_time_scaled=shot_time/range(shot_time);
cdat=[interp1(c_cord,cdat(:,1),shot_time_scaled),...
    interp1(c_cord,cdat(:,2),shot_time_scaled),...
    interp1(c_cord,cdat(:,3),shot_time_scaled)];

scatter((probe_freq-nanmean(probe_freq))*1e-9,corrected_delta,30,cdat,'square','filled')
colormap(viridis(1000))
c =colorbar;
c.Label.String = 'time (H)';
caxis([0,range(shot_time)/(60*60)])
xlabel('delta probe beam frequency (GHz)')
ylabel('delta freq')
subplot(2,1,2)
square_diff= (3*(1/anal_opts.atom_laser.pulsedt)+...
    data.osc_fit.model_coefs(probe_dat_mask,2,1)).^2-...
    (freq_drift_model(data.mcp_tdc.time_create_write(probe_dat_mask,1))).^2;
scatter((probe_freq-nanmean(probe_freq))*1e-9,square_diff,30,cdat,'square','filled')
xlabel('delta probe beam frequency (GHz)')
ylabel('square difference in freq')
colormap(viridis(1000))
c =colorbar;
c.Label.String = 'time (H)';
caxis([0,range(shot_time)/(60*60)])

set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1600, 900])
data.mcp_tdc.all_ok=tmp_all_ok;
plot_name='tuneout_time_graph';
saveas(gcf,[anal_out.dir,plot_name,'.png'])
saveas(gcf,[anal_out.dir,plot_name,'.fig'])



%%

%%
%now we do some fitting
%thresholds for CI
%sd         CI
%1          0.3174
%2          0.05
%3          2.699e-03
ci_size_disp=0.3174;%one sd %confidence interval to display
ci_size_cut_outliers=0.01; %confidence interval for cutting outliers

fprintf('Calculating Fits')
%select the data in some freq range and that has an ok number

%set up the data input for the fit
xdat=probe_freq(~isnan(probe_freq))';
ydat=square_diff(~isnan(probe_freq))';
cdat=cdat(~isnan(probe_freq));
if exist('new_to_freq_val','var') %if the TO has been calculated before use that as the center
    freq_offset=new_to_freq_val;
else
    freq_offset=nanmean(xdat); %otherwise use the center of the range
end
xdat=xdat-freq_offset; %fits work better when thery are scaled reasonably
xdat=xdat*1e-9;

modelfun = @(b,x) b(1)+ b(2).*x; %simple linear model
opts = statset('nlinfit');
%opts.RobustWgtFun = 'welsch' ; %a bit of robust fitting
%opts.Tune = 1;
beta0 = [1e-5,1e-2]; %intial guesses
mdl_all = fitnlm(xdat,ydat,modelfun,beta0,'Options',opts);
xsamp=linspace(min(xdat),max(xdat),1e3)'; %sample for the model curve
[ysamp,yci]=predict(mdl_all,xsamp,'Prediction','observation','Alpha',ci_size_disp); %note the observation CI
%now plot the data and the model together
figure(6);
set(gcf,'color','w')
subplot(1,2,1)
plot(xsamp,ysamp,'k-')
hold on
plot(xsamp,yci,'r-')
scatter(xdat,ydat,30,cdat,'square','filled')
c =colorbar;
c.Label.String = 'time (H)';
caxis([0,range(shot_time)/(60*60)])
%plot(xdat,ydat,'bx')
xlabel(sprintf('probe beam set freq - %.3f(GHz)',freq_offset*1e-9))
ylabel('Response (Hz^2)')
title('Good Data')
first_plot_lims=[get(gca,'xlim');get(gca,'ylim')];
%now make a new prediction with the model but with the CI to cut out outliers
[~,yci_cull_lim]=predict(mdl_all,xdat','Prediction','observation','Alpha',ci_size_cut_outliers);
is_outlier_idx=ydat>yci_cull_lim(:,1)' & ydat<yci_cull_lim(:,2)';
%color the ones that will be removed
plot(xdat(~is_outlier_idx),ydat(~is_outlier_idx),'r.','markersize',15)
hold off
xdat_culled=xdat(is_outlier_idx);
ydat_culled=ydat(is_outlier_idx);
cdat_culled=cdat(is_outlier_idx);

opts = statset('nlinfit');
%opts.RobustWgtFun = 'welsch' ;
opts.Tune = 1;
beta0 = [1e-5,1e-2];
mdl_culled = fitnlm(xdat_culled,ydat_culled,modelfun,beta0,'Options',opts);
xsamp_culled=linspace(min(xdat_culled),max(xdat_culled),1e3)';
[ysamp_culled,yci_culled]=predict(mdl_culled,xsamp_culled,'Alpha',0.2); %'Prediction','observation'
%now plot the remaining data along with the fit model and the model CI


subplot(1,2,2)
plot(xsamp_culled,ysamp_culled,'k-')
hold on
plot(xsamp_culled,yci_culled,'r-')
scatter(xdat_culled,ydat_culled,30,cdat_culled,'square','filled')
colormap(viridis(1000))
c =colorbar;
c.Label.String = 'time (H)';
caxis([0,range(shot_time)/(60*60)])
hold off
xlabel(sprintf('probe beam set freq - %.3f (GHz)',freq_offset*1e-9))
ylabel('Response (Hz^2)')
title('Fit Outliers Removed')
set(gca,'xlim',first_plot_lims(1,:))
set(gca,'ylim',first_plot_lims(2,:))


set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1600, 900])
data.mcp_tdc.all_ok=tmp_all_ok;
plot_name='TO_fits';
saveas(gcf,[anal_out.dir,plot_name,'.png'])
saveas(gcf,[anal_out.dir,plot_name,'.fig'])



%normalize by the CI at the TO
figure(7);
clf
set(gcf,'color','w')
cross_xval=-mdl_all.Coefficients.Estimate(1)/mdl_all.Coefficients.Estimate(2);
[cross_yval,cross_yci]=predict(mdl_all,cross_xval,'Prediction','observation','Alpha',ci_size_disp);
if abs(cross_yval)>1e-2, error('not crossing zero here') ,end
cross_yci=diff(cross_yci)/2;
plot(xsamp,ysamp/cross_yci,'k-')
hold on
plot(xsamp,yci/cross_yci,'r-')
plot(xdat,ydat/cross_yci,'bx')
hold off
xlabel(sprintf('probe beam set freq - %.3f (GHz)',freq_offset*1e-9))
ylabel('Response scaled to sample SD')
title('Senistivity Graph ')

set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1600, 900])
data.mcp_tdc.all_ok=tmp_all_ok;
plot_name='Sens_graph';
saveas(gcf,[anal_out.dir,plot_name,'.png'])
saveas(gcf,[anal_out.dir,plot_name,'.fig'])


%inverse scaled gradient to give the single shot uncert
single_shot_uncert=abs(1/(mdl_all.Coefficients.Estimate(2)/cross_yci));
fprintf('\nfit results\n')
fprintf('median damping time %.2f\n',median(1./data.osc_fit.model_coefs(data.osc_fit.ok.rmse,7,1)))
%calculate some statistics and convert the model parameter into zero crossing and error therin
old_to_wav=413.0938e-9;
new_to_freq_val=-1e9*mdl_culled.Coefficients.Estimate(1)/mdl_culled.Coefficients.Estimate(2)+freq_offset;
new_to_freq_unc=1e9*abs((mdl_culled.Coefficients.Estimate(1)/mdl_culled.Coefficients.Estimate(2)))*...
    sqrt((mdl_culled.Coefficients.SE(1)/mdl_culled.Coefficients.Estimate(1))^2+...
    (mdl_culled.Coefficients.SE(2)/mdl_culled.Coefficients.Estimate(2))^2);
new_to_wav_val=const.c/(new_to_freq_val*2);
new_to_wav_unc=new_to_freq_unc*const.c/((new_to_freq_val*2)^2);
fprintf('run start time %.1f\n (posix)',...
    data.mcp_tdc.time_create_write(1,2)-anal_opts.trig_dld-anal_opts.dld_aquire)
fprintf('run stop time  %.1f\n (posix)',...
    data.mcp_tdc.time_create_write(end,2)-anal_opts.trig_dld-anal_opts.dld_aquire)
fprintf('duration  %.1f\n (s)',...
    data.mcp_tdc.time_create_write(end,2)-data.mcp_tdc.time_create_write(1,2))
fprintf('TO freq        %.2f±%.2f MHz\n',new_to_freq_val*1e-6,new_to_freq_unc*1e-6)
fprintf('TO wavelength  %.6f±%f nm \n',new_to_wav_val*1e9,new_to_wav_unc*1e9)
fprintf('diff from TOV1 %e±%e nm \n',(new_to_wav_val-old_to_wav)*1e9,new_to_wav_unc*1e9)
fprintf('number of probe files        %u \n',sum(probe_dat_mask))
fprintf('number of calibration files  %u \n',sum(cal_dat_mask))
fprintf('total used                   %u \n',sum(probe_dat_mask)+sum(cal_dat_mask))
fprintf('files with enough number     %u\n',sum(data.mcp_tdc.num_ok))

fprintf('single shot uncert detuning @1SD %.1f MHz, %.2f fm\n',single_shot_uncert*1e3,...
    1e9*single_shot_uncert*const.c/((new_to_freq_val*2)^2)*10^15)
%predicted uncert using this /sqrt(n)
fprintf('predicted stat. uncert %.1f MHz, %.2f fm\n',single_shot_uncert/sqrt(sum(data.mcp_tdc.num_ok))*1e3,...
    1e9*single_shot_uncert/sqrt(sum(data.mcp_tdc.num_ok))*const.c/((new_to_freq_val*2)^2)*10^15)

diary off

%% try and find what the outliers were doing
%given this mask ~color_idx find the shot nums and times 
%~isnan(probe_freq)


%% damping results
%plot out what the distibution over damping times is
% figure(7);
% set(gcf,'color','w')
% histogram(1./data.osc_fit.model_coefs(data.osc_fit.ok.rmse,7,1),linspace(0,3,1e2))

toc