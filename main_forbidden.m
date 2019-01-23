% measure the 427nm forbidden transtion


%monitor the forbidden transtion run
% monitor an experiment that uses measrument of the trap freq
% It would be usefull to get a decent idea of what the trap frequency is doing during a run without
% the need for a full processing of the data as in main
%  - processing each shot as it is made
%  - plot a history of the trap freq

% Known BUGS/ Possible Improvements
%
%
% Author: Bryce Henson
% email: Bryce.Henson@live.com
% Last revision:2018-10-01
% BEGIN USER VAR-------------------------------------------------
anal_opts.tdc_import.dir='Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20190102_427nm_run_prelim\';
anal_opts.tdc_import.file_name='d';
anal_opts.tdc_import.force_load_save=false;   %takes precidence over force_reimport
anal_opts.tdc_import.force_reimport=true;
anal_opts.tdc_import.force_forc=false;
anal_opts.tdc_import.dld_xy_rot=0.61;

tmp_xlim=[-30e-3, 30e-3];     %tight XY lims to eliminate hot spot from destroying pulse widths
tmp_ylim=[-30e-3, 30e-3];
tlim=[0,4];
anal_opts.tdc_import.txylim=[tlim;tmp_xlim;tmp_ylim];

anal_opts.max_runtime=inf;%cut off the data run after some number of hours
anal_opts.atom_laser.pulsedt=8.000e-3;
anal_opts.atom_laser.t0=0.41784; %center i ntime of the first pulse
anal_opts.atom_laser.start_pulse=1; %atom laser pulse to start with
anal_opts.atom_laser.pulses=325;
anal_opts.atom_laser.appr_osc_freq_guess=[52,48,40];
anal_opts.atom_laser.pulse_twindow=anal_opts.atom_laser.pulsedt*0.9;
anal_opts.atom_laser.xylim=anal_opts.tdc_import.txylim(2:3,:); %set same lims for pulses as import

anal_opts.global.fall_time=0.417;
anal_opts.global.qe=0.09;

anal_opts.trig_dld=20.3;
anal_opts.dld_aquire=4;
anal_opts.trig_ai_in=20;
anal_opts.ai_log.probe_set_pt=1.5;



% END USER VAR-----------------------------------------------------------
%sets up the struct 'data' which will contain everything you could want incuding the txy data and
%the information from the logs
data=[]; %CLEAR THE DATA
anal_out=[];
%set up an output dir %https://gist.github.com/ferryzhou/2269380
if anal_opts.tdc_import.dir(end) ~= filesep, anal_opts.tdc_import.dir = [anal_opts.tdc_import.dir filesep]; end
if (exist([anal_opts.tdc_import.dir,'out'], 'dir') == 0), mkdir([anal_opts.tdc_import.dir,'out']); end
 
anal_out.dir=sprintf('%sout\\%s\\',...
    anal_opts.tdc_import.dir,datestr(datetime('now'),'yyyymmddTHHMMSS'));
if (exist(anal_out.dir, 'dir') == 0), mkdir(anal_out.dir); end
anal_opts.global.out_dir=anal_out.dir;
 
diary([anal_out.dir,'anal.txt'])

%add all subfolders to the path
this_folder = fileparts(which(mfilename));
% Add that folder plus all subfolders to the path.
addpath(genpath(this_folder));

hebec_constants %call the constants function that makes some globals
anal_opts.global.velocity=const.g0*anal_opts.global.fall_time;

%%

%% IMPORT TDC DATA to data.mcp_tdc
anal_opts.tdc_import.shot_num=find_data_files(anal_opts.tdc_import);
%mask= anal_opts.tdc_import.shot_num<=1700;
%anal_opts.tdc_import.shot_num=anal_opts.tdc_import.shot_num(mask);
%anal_opts.tdc_import.shot_num= anal_opts.tdc_import.shot_num(1:10); %debuging
[mcp_tdc_data,import_opts]=import_mcp_tdc_data(anal_opts.tdc_import);
data.mcp_tdc=mcp_tdc_data;
%% IMPORT LV LOG to data.labview
%TO DO FUNCTIONALIZE
%import the wavemeter log
%adaptively to deal with the 2 different log files that are in the data
lv_log=[];
lv_log.dir = strcat(anal_opts.tdc_import.dir,'log_LabviewMatlab.txt');
fid = fopen(lv_log.dir );
lv_log.cell=textscan(fid,'%s','Delimiter','\n');
fclose(fid);
lv_log.cell=lv_log.cell{1};
for ii=1:size(lv_log.cell,1)
    if ~isequal(lv_log.cell{ii},'') %catch the empty case
        if contains(lv_log.cell{ii},'measure_probe')
            line_cells=textscan(lv_log.cell{ii},'%f %s %s %s %f %s %u','Delimiter',',');
            lv_log.setpoints(ii)=line_cells{5};
            lv_log.probe_calibration(ii)=false;
            lv_log.iter_nums(ii)=line_cells{7};
        elseif contains(lv_log.cell{ii},'calibrate')
            line_cells=textscan(lv_log.cell{ii},'%f %s %s %s %s %u','Delimiter',',');
            lv_log.setpoints(ii)=NaN;
            lv_log.probe_calibration(ii)=true;
            lv_log.iter_nums(ii)=line_cells{6};
        else %deals with the legacy case (only 20180813_CW_AL_tuneout_scan)
            line_cells=textscan(lv_log.cell{ii},'%f %s %s %s %f %s %u','Delimiter',',');
            lv_log.setpoints(ii)=line_cells{5};
            lv_log.probe_calibration(ii)=false;
            lv_log.iter_nums(ii)=line_cells{7};
        end
        lv_log.posix_times(ii)=line_cells{1};
        lv_log.iso_times{ii}=line_cells{2};
    end
end
data.labview=[];
data.labview.setpoint=lv_log.setpoints*1e6; %convert to hz
data.labview.time=lv_log.posix_times;
data.labview.shot_num=lv_log.iter_nums;
data.labview.calibration=lv_log.probe_calibration;
%can check that the times look ok
% plot(data.mcp_tdc.write_time-data.probe.time)

%% CHECK ATOM NUMBER
%total number of detected counts
sfigure(1);
clf
set(gcf,'color','w')
subplot(4,1,1)
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
%should plot the threshold

%% Match Labview data
%because the mcp-dld detector has no direct communication with the bec computer
% the data.labview.shot_num does not nessesarily correspond to data.mcp_tdc.shot_num
%try and match up the file with if it is a calibaration using the time
%it is slightly overkill here to search each one, but just being extra
%cautious/flexible
time_thresh=4; %how close for the times to be considered the same shot
%lets examine what the time difference does

data.mcp_tdc.labview_shot_num=[];
data.mcp_tdc.probe.calibration=[];

imax=min([size(data.labview.time,2),size(data.mcp_tdc.time_create_write,1)]);
%imax=5000;
time_diff=data.mcp_tdc.time_create_write(1:imax,2)'-anal_opts.dld_aquire-anal_opts.trig_dld-...
    data.labview.time(1:imax);
mean_delay_labview_tdc=median(time_diff);

sfigure(1);
set(gcf,'color','w')
subplot(4,1,2)
plot(data.mcp_tdc.shot_num,time_diff)
xlabel('shot number')
ylabel('time between labview and mcp tdc')
title('raw time diff')
%to do include ai_log
iimax=size(data.mcp_tdc.time_create_write(:,1),1);
data.mcp_tdc.probe.calibration=nan(iimax,1);
data.mcp_tdc.labview_shot_num=nan(iimax,1);
%loop over all the tdc_files 
for ii=1:iimax
    %predict the labview master trig time
    %use the write time to handle being unpacked from 7z
    est_labview_start=data.mcp_tdc.time_create_write(ii,2)...
        -anal_opts.trig_dld-anal_opts.dld_aquire-mean_delay_labview_tdc;
    [tval,nearest_idx]=closest_value(data.labview.time...
        ,est_labview_start);
    if abs(tval-est_labview_start)<time_thresh
        data.mcp_tdc.labview_shot_num(ii)=data.labview.shot_num(nearest_idx);
        data.mcp_tdc.probe.calibration(ii)=data.labview.calibration(nearest_idx);
    end 
end


%% IMPORT THE ANALOG INPUT LOG
%load in the analog input files to check if the laser is single mode & that the potodiode value is close to the set
%point
% a two teired cache system is used one level for importing all 

anal_opts.ai_log.dir=anal_opts.tdc_import.dir;
anal_opts.ai_log.force_reimport=false;
anal_opts.ai_log.force_load_save=false;
anal_opts.ai_log.log_name='log_analog_in_';
anal_opts.ai_log.pd.set=data.mcp_tdc.probe.calibration;
%nan compatable logical inverse
anal_opts.ai_log.pd.set(~isnan(anal_opts.ai_log.pd.set))=~anal_opts.ai_log.pd.set(~isnan(anal_opts.ai_log.pd.set));
anal_opts.ai_log.pd.set=anal_opts.ai_log.pd.set*anal_opts.ai_log.probe_set_pt;
%anal_opts.ai_log.pd.set(isnan(anal_opts.ai_log.pd.set))=0;
anal_opts.ai_log.aquire_time=4;

anal_opts.ai_log.pd.diff_thresh=0.1;
anal_opts.ai_log.pd.std_thresh=0.1;
anal_opts.ai_log.pd.time_start=0.44;
anal_opts.ai_log.pd.time_stop=1.225;
anal_opts.ai_log.sfp.num_checks=10; %how many places to check that the laser is single mode
anal_opts.ai_log.sfp.thresh_cmp_peak=20e-3; %theshold on the compressed signal to be considered a peak
anal_opts.ai_log.sfp.peak_dist_min_pass=4.5;%minimum (min difference)between peaks for the laser to be considered single mode
anal_opts.ai_log.plot.all=false;
anal_opts.ai_log.plot.failed=false;
anal_opts.ai_log.time_match_valid=5; %how close the predicted start of the shot is to the actual
anal_opts.ai_log.scan_time=14e-3;  %estimate of the sfp scan time,used to set the window and the smoothing


%because im only passing the ai_log feild to aviod conflicts forcing a reimport i need to coppy these feilds
anal_opts.ai_log.trig_dld=anal_opts.trig_dld;
anal_opts.ai_log.dld_aquire=anal_opts.dld_aquire;
anal_opts.ai_log.aquire_time=anal_opts.dld_aquire;
anal_opts.ai_log.trig_ai_in=anal_opts.trig_ai_in;
ai_log_out=ai_log_import(anal_opts.ai_log,data);
%copy the output across
data.ai_log=ai_log_out;

% % Trying to automate setpoint correction

%save([datestr(datetime('now'),'yyyymmddTHHMMSS'),'.mat'],'-v7.3')

%% IMPORT WM LOG FILES

anal_opts.wm_log.dir=anal_opts.tdc_import.dir;

anal_opts.wm_log.force_reimport=false;
wm_log_name='log_wm_';
wm_logs=dir([anal_opts.wm_log.dir,wm_log_name,'*.txt']);
anal_opts.wm_log.names={wm_logs.name};
data.wm_log.raw=wm_log_import(anal_opts.wm_log);


%% CHECK THE WM INPUTS
%check that the probe beam (optical) freq was stable & that 2x red ~ blue
%(now redundant) check that the doubler photodiode voltage is ok
%define a time window for checking if the doubler was ok &averaging the wavelength of the laser
%i think it will be anal_opts.atom_laser.t0 after the creation time of the tdc file
%compexity is that the time that the tdc file is wrote/reated is not relaible and depend on the flux rate and avaialble mem
%to this end find the closest labview update time and go back one then fowards

anal_opts.wm_log.plot_all=false;
anal_opts.wm_log.plot_failed=false;
anal_opts.wm_log.force_reimport=false;

anal_opts.wm_log.time_pd_padding=4; %check this many s each side of probe
anal_opts.wm_log.time_blue_padding=1; %check this many seconde each side of probe
anal_opts.wm_log.time_probe=3;
anal_opts.wm_log.ecd_volt_thresh=0.5;

anal_opts.wm_log.red_sd_thresh=1; %allowable standard deviation in MHz
anal_opts.wm_log.red_range_thresh=1; %allowable range deviation in MHz
anal_opts.wm_log.rvb_thresh=10; %allowable value of abs(2*red-blue)

data.wm_log.proc=wm_log_process(anal_opts,data);
clear('sub_data')

%%
%% COMBINE ALL CHECK LOGICS AND PLOT
%Here we will do a plot of all the checks and then combine them into one
%master 'ok'/check vector
%here we keep this vector of ok logic the same size as the data.mcp_tdc to simplify use later

%data.mcp_tdc.probe.ok.reg_pd
%data.mcp_tdc.probe.ok.sfp
%data.mcp_tdc.num_ok
%data.mcp_tdc.probe.ok.freq; %frequency reading is stable
%data.mcp_tdc.probe.ok.rvb;  %2r-b check
%data.mcp_tdc.probe.ok.ecd_pd;  %ecd pd value
  

sfigure(1);
set(gcf,'color','w')
subplot(4,1,3)
%plot all the logics, dither it a bit to make it easier to figure out
%culprits
line_width=1.5;
stairs(data.mcp_tdc.shot_num,data.mcp_tdc.num_ok-0.03,'LineWidth',line_width)
hold on


stairs(data.mcp_tdc.shot_num,data.ai_log.ok.reg_pd-0.01,'LineWidth',line_width)
stairs(data.mcp_tdc.shot_num,data.ai_log.ok.sfp-0.02,'LineWidth',line_width)
stairs(data.mcp_tdc.shot_num,data.wm_log.proc.ok.freq+0.00,'LineWidth',line_width)
%forms another check that the laser is single mode
stairs(data.mcp_tdc.shot_num,data.wm_log.proc.ok.rvb+0.01,'LineWidth',line_width)
 %this is reduncant as it should be caught by the reg_pd measurement
stairs(data.mcp_tdc.shot_num, data.wm_log.proc.ok.ecd_pd+0.02,'LineWidth',line_width) 
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

subplot(4,1,4)
tmp_cal=data.mcp_tdc.probe.calibration;
tmp_cal(isnan(tmp_cal))=false;
%must have good atom number AND (good probe OR be calibration)
tmp_probe_ok=(data.ai_log.ok.sfp &...
    data.ai_log.ok.reg_pd &...
    data.wm_log.proc.ok.freq &....
    data.wm_log.proc.ok.rvb );
tmp_all_ok=data.mcp_tdc.num_ok' &...
    (tmp_probe_ok| tmp_cal);
    %data.mcp_tdc.probe.ok.ecd_pd;
stairs(data.mcp_tdc.shot_num,tmp_probe_ok+0.01,'LineWidth',line_width)
hold on
stairs(data.mcp_tdc.shot_num,tmp_all_ok,'LineWidth',line_width)
hold off
legend('all','probe')
ylabel('Good?')
xlabel('Shot Number')
set(gca,'ytick',[0,1],'yticklabel',{'False','True'})
title('ALL ok')
ylim([-0.1,1.1])
tmp_num_shots=numel(data.mcp_tdc.shot_num);
tmp_num_ok_shots=sum(tmp_all_ok);
data.mcp_tdc.all_ok=tmp_all_ok;
data.mcp_tdc.probe.ok.all=tmp_probe_ok;


fprintf('ok logic gives %u / %u shots for yeild %04.1f %%\n',...
    tmp_num_ok_shots,tmp_num_shots,1e2*tmp_num_ok_shots/tmp_num_shots)
% set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1600, 900])
% plot_name='check_logics';
% saveas(gcf,[anal_out.dir,plot_name,'.png'])
% saveas(gcf,[anal_out.dir,plot_name,'.fig'])

%%

%%
%% BINNING UP THE ATOM LASER PULSES
%now find the mean position of each pulse of the atom laser in each shot
data.mcp_tdc.al_pulses=bin_al_pulses(anal_opts.atom_laser,data);

%% Summing the AOM number
pulses_seg1=[1,15];
pulses_seg2=[144,210];

tot_counts_seg1=nansum(data.mcp_tdc.al_pulses.num_counts(:,pulses_seg1(1):pulses_seg1(2)),2);
tot_counts_seg2=nansum(data.mcp_tdc.al_pulses.num_counts(:,pulses_seg2(1):pulses_seg2(2)),2);
mask=data.mcp_tdc.all_ok;

sfigure(3);
set(gcf,'color','w')
plot(probe_freqs(mask),...
    tot_counts_seg1(mask)./tot_counts_seg2(mask),...
    'k.')
xlabel('Shot Number')
ylabel('Log Number remaining')
pause(1e-6)



%% FITTING THE ATOM NUMBER PART 1
%use the inital few atom laser pulses in order to determine the atom number
%not of that much benifit TBH
anal_opts.atom_num_fit=[];
anal_opts.atom_num_fit.pulses=pulses_seg1; %min,max index of pulses
sfigure(2)
subplot(2,1,1)
%only show the fit if there is less than 5 shots
anal_opts.atom_num_fit.plot.each_shot=false;
anal_opts.atom_num_fit.plot.history=false;
anal_opts.atom_num_fit.qe=anal_opts.global.qe;

data.num_fit.pre_probe=fit_atom_number(anal_opts.atom_num_fit,data);

%% FITTING THE ATOM NUMBER PART 2
sfigure(2)
subplot(2,1,2)
anal_opts.atom_num_fit.plot.each_shot=false;
anal_opts.atom_num_fit.pulses=pulses_seg2; %min,max index of pulses
data.num_fit.post_probe=fit_atom_number(anal_opts.atom_num_fit,data);

%% Compare
mask=data.mcp_tdc.all_ok;

atom_num_pre=cellfun(@(x) x(pulses_seg1(2)+1),data.num_fit.pre_probe.fit_predict(mask));
atom_num_post=cellfun(@(x) x(pulses_seg2(1)),data.num_fit.post_probe.fit_predict(mask));
atom_num_unc_pre=cellfun(@(x) x(pulses_seg1(2)+1),data.num_fit.pre_probe.fit_predict_unc(mask));
atom_num_unc_post=cellfun(@(x) x(pulses_seg2(1)),data.num_fit.post_probe.fit_predict_unc(mask));

frac_rem_vec_lin=atom_num_post./atom_num_pre;
frac_rem_vec_log_unc=frac_rem_vec_lin.*sqrt((atom_num_unc_pre./atom_num_pre).^2+(atom_num_unc_post./atom_num_post).^2);
frac_rem_vec_log_val=log(atom_num_post./atom_num_pre);
frac_rem_vec_log_unc=abs(sqrt((atom_num_unc_pre./atom_num_pre).^2+(atom_num_unc_post./atom_num_post).^2));
probe_freqs=data.wm_log.proc.probe.freq.act.mean(mask);

%%
sfigure(3);
set(gcf,'color','w')
errorbar(probe_freqs,...
    frac_rem_vec_lin,frac_rem_vec_log_unc,...
    'k.','MarkerSize',7,'CapSize',0,'LineWidth',1)
xlabel('Shot Number')
ylabel('Log Number remaining')
pause(1e-6)


%%
sfigure(4);
set(gcf,'color','w')
plot(data.wm_log.proc.probe.freq.act.mean,data.mcp_tdc.num_counts,'k.')
