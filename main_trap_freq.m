%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine the Tune out from a dataset
% using the measured change in trap frequency from the application of a probe beam.
% application of the tune out probe beam.
% The script:
% - defines the user controled options
% - Imports all the tdc data files 
% - Imports labview log file
% - Imports the wavemeter log file
% - Match up the Labview data withthe tdc data
% - Imports the analog in log file , for each file the import:
%   - checks that the pd voltage is ok
%   - check that the laser is single mode using the scanning fabry perot signals
%   - fits the measured AC mains waveform
% - Check that the wavemeter readings are ok for each shot
%   - checks that the wavelengths is stable during the probe intterogation
%   - checks that the red wavelength is ~half the blue
%   - checks that the double pd voltage is ok (now redundant beacuse of probe pd)
% - checks that the number of counts in the file is ok
% - combines all these checks into one master check
% - bins up each pulse of the AL
% - Fits the trap frequency
%   - Investigate fit correlations
% - mask out only the (good)calibrations shots and make a model of how the (unpeturbed) trap freq changes in time
% - calculate the probe beam signal
% - fit the relation of the probe beam signal to the probe optical frequency
%   - determine the TO
%   - bootstrap the output uncert 
%   - fit with a linear and higher order model
% - fit the relation of the probe bem signal to the probe optical frequency for each scan
%   - investigate correlations

%the data structure
%   first level is instrument or anal method
%   try and pass things between the modules of code only using the 'data' struct
%   keep evertthing referenced to data.mcp_tdc eg. the probe beam power ok should corespond to 

% TIMING 
% because there are so many moving peices a lot of the script requires matching up the times of the varrious inputs
% LABVIEW writes to the log
%    | (~0.25s)
% DAC master trig ---------->		Digital output cards
%									|			|
%									|			|(anal_opts.trig_ai_in ~20s)
%	(anal_opts.trig_dld~20.3s)		|			|
%									|			Trig analog in
%					DLD trig,dld create file	|
%									|			|(analog in aq time)
%			(dld aq time)			|			analog in end
%									|			file creation,file modified
%								DLD write
%

% Other m-files required: import_data,find_data_files,dld_raw_to_txy,masktxy,data_tcreate,
%                         dld_read_5channels_reconst_multi_imp,constants
% Also See:
% Subfunctions: none
% MAT-files required: none
%
% Known BUGS/ Possible Improvements
%  - a better way of figuring out what the pd setpoint is, maybe write as a config file
%  - harmonize the output
%   -the fit error depends on wavelength indicating that the model does not
%   have enough freedom
%	-make plots more compact
%   -harmonize the anal opts

%       
%
% Author: Bryce Henson
% email: Bryce.Henson@live.com
% Last revision:2019-05-26



%close all
clear all
%%
% BEGIN USER VAR-------------------------------------------------

%setup directories you wish to loop over
%% for testing
% loop_config.dir = {
%     '..\scratch_data\20190227_qwp_270',
%     '..\scratch_data\20190227_qwp_286',
%     '..\scratch_data\20190227_qwp_310',g
%     };
% loop_config.set_pt = [nan,nan,nan];
% selected_dirs = 1:numel(loop_config.dir); %which files to loop over (currently all)
% selected_dirs = [1,2,3];

%% for deployment
% select the directories in a folder
root_data_dir='G:\good_data';
%root_data_dir='..\scratch_data';
files = dir(root_data_dir);
files=files(3:end);
% Get a logical vector that tells which is a directory.
dir_mask = [files.isdir];
folders=files(dir_mask);
%convert to the full path
folders=arrayfun(@(x) fullfile(root_data_dir,x.name),folders,'UniformOutput' ,false);
loop_config.dir = folders;
loop_config.set_pt=nan(1,numel(loop_config.dir));
selected_dirs=1:numel(loop_config.dir);

%%
anal_opts=[]; %reset the options (would be good to clear all variables except the loop config
anal_opts.tdc_import.file_name='d';
anal_opts.tdc_import.force_load_save=false;   %takes precidence over force_reimport
anal_opts.tdc_import.force_reimport=false;
anal_opts.tdc_import.force_forc=false;
anal_opts.tdc_import.dld_xy_rot=0.61;
anal_opts.tdc_import.save_cache_in_data_dir=true;

%Should probably try optimizing these
tmp_xlim=[-30e-3, 30e-3];     %tight XY lims to eliminate hot spot from destroying pulse widths
tmp_ylim=[-30e-3, 30e-3];
tlim=[0,4];
anal_opts.tdc_import.txylim=[tlim;tmp_xlim;tmp_ylim];
anal_opts.max_runtime=inf;%inf%cut off the data run after some number of hours, should bin included as its own logic not applied to the atom number ok

anal_opts.trig_dld=20.3;
anal_opts.dld_aquire=4;
anal_opts.aquire_time=4;
anal_opts.trig_ai_in=20;
anal_opts.aom_freq= 189.*1e6;%Hz %set to zero for comparison with previous data runs

anal_opts.wm_log.plot_all=false;
anal_opts.wm_log.plot_failed=false;


anal_opts.atom_laser.t0=0.417770; %center i ntime of the first pulse

anal_opts.global.fall_time=0.417;
anal_opts.global.qe=0.09;
anal_opts.global.atom_laser.t0=anal_opts.atom_laser.t0;


%anal_opts.osc_fit.binsx=1000;
%anal_opts.osc_fit.blur=1;
%anal_opts.osc_fit.xlim=[-20,20]*1e-3;
%anal_opts.osc_fit.tlim=[0.86,1.08];

date_str='20190601T000000';
reprocess_folder_if_older_than=posixtime(datetime(datenum(date_str,'yyyymmddTHHMMSS'),'TimeZone','local','ConvertFrom','datenum'));%posix date
active_process_mod_time=60*5;

% END USER VAR-----------------------------------------------------------


%% set up the path

% find this .m file's path, this must be in the project root dir
this_folder = fileparts(which(mfilename));
% Add that folder plus all subfolders to the path.
addpath(genpath(this_folder));%add all subfolders to the path to find genpath_exclude
path_to_genpath=fileparts(which('genpath_exclude'));
path(pathdef) %clean up the path back to the default state to remove all the .git that were added
addpath(this_folder)
addpath(path_to_genpath)
addpath(genpath_exclude(fullfile(this_folder,'lib'),'\.')) %dont add hidden folders
addpath(genpath_exclude(fullfile(this_folder,'dev'),'\.'))
addpath(genpath_exclude(fullfile(this_folder,'bin'),'\.'))

hebec_constants %call the constants function that makes some globals

% loop over the selected directories
for dir_idx = selected_dirs
anal_opts.tdc_import.dir = loop_config.dir{dir_idx};
%check that this folder is not being processed by another matlab instance
if should_process_folder(anal_opts.tdc_import.dir,reprocess_folder_if_older_than,active_process_mod_time)  
try
    
fprintf('processing data dir %s \n',anal_opts.tdc_import.dir)
main_trap_freq_timer=tic;
anal_opts.probe_set_pt=loop_config.set_pt(dir_idx);
%sets up the struct 'data' which will contain everything you could want incuding the txy data and
%the information from the logs
data=[]; %CLEAR THE DATA
anal_out=[];
%add a file seperator to the end of the import path
if anal_opts.tdc_import.dir(end) ~= filesep, anal_opts.tdc_import.dir = [anal_opts.tdc_import.dir filesep]; end

anal_opts.global.fall_velocity=const.g0*anal_opts.global.fall_time; %velocity when the atoms hit the detector
% fall_dist=1/2 a t^2 
%TODO get from engineering documents
anal_opts.global.fall_dist=(1/2)*const.g0*anal_opts.global.fall_time^2;

%% IMPORT TDC DATA to data.mcp_tdc
anal_opts.tdc_import.shot_num=find_data_files(anal_opts.tdc_import);
%anal_opts.tdc_import.shot_num= anal_opts.tdc_import.shot_num(1:10); %debuging

%set up an output dir %https://gist.github.com/ferryzhou/2269380
if (exist([anal_opts.tdc_import.dir,'out'], 'dir') == 0), mkdir([anal_opts.tdc_import.dir,'out']); end
%make a subfolder with the ISO timestamp for that date
anal_out.dir=sprintf('%sout\\%s\\',...
    anal_opts.tdc_import.dir,datestr(datetime('now'),'yyyymmddTHHMMSS'));
if (exist(anal_out.dir, 'dir') == 0), mkdir(anal_out.dir); end
anal_opts.global.out_dir=anal_out.dir;
%start up the diary of stdout
diary([anal_out.dir,'anal.txt'])
%import the data
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

%create a list of indicies (of the mcp_tdc) that have an ok number of counts
%exclude the very low and then set the thresh based on the sd of the remaining
not_zero_files=data.mcp_tdc.num_counts>1e3; 
num_thresh=0.5*median(data.mcp_tdc.num_counts(not_zero_files));
data.mcp_tdc.num_ok=data.mcp_tdc.num_counts>num_thresh & ...
    (data.mcp_tdc.time_create_write(:,1)'-data.mcp_tdc.time_create_write(1,1))<(anal_opts.max_runtime*60*60);
fprintf('shots number ok %u out of %u \n',sum(data.mcp_tdc.num_ok),numel(data.mcp_tdc.num_ok))
drawnow
%plot((data.mcp_tdc.time_create_write(:,2)-data.mcp_tdc.time_create_write(1,2))/(60*60),data.mcp_tdc.num_counts)
%xlabel('time (h)')
stfig('diagnostics and veto','add_stack',1);
clf
subplot(4,1,1)
plot(data.mcp_tdc.shot_num,data.mcp_tdc.num_counts,'k')
xlabel('shot number')
ylabel('total counts')
title('num count run trend')
%should plot the threshold

%% Match Labview data
%because the mcp-dld detector has no direct communication with the bec computer
% the data.labview.shot_num does not nessesarily correspond to data.mcp_tdc.shot_num
%try and match up the file with if it is a calibaration using the time
%it is slightly overkill here to search each one, but just being extra
%cautious/flexible
time_thresh=3; %4how close for the times to be considered the same shot
%lets examine what the time difference does

data.mcp_tdc.labview_shot_num=[];
data.mcp_tdc.probe.calibration=[];

imax=min([size(data.labview.time,2),size(data.mcp_tdc.time_create_write,1)]);
%imax=5000;
time_diff=data.mcp_tdc.time_create_write(1:imax,2)'-anal_opts.dld_aquire-anal_opts.trig_dld-...
    data.labview.time(1:imax);
mean_delay_labview_tdc=0;%median(time_diff);

stfig('diagnostics and veto','add_stack',1);
subplot(4,1,2)
plot(data.mcp_tdc.shot_num(1:imax),time_diff-mean_delay_labview_tdc,'k')
xlabel('shot number')
ylabel('corrected time between labview and mcp tdc')
title('raw time diff')
drawnow
%to do include ai_log
iimax=size(data.mcp_tdc.time_create_write(:,1),1);
data.mcp_tdc.probe.calibration=nan(iimax,1);
data.mcp_tdc.labview_shot_num=nan(iimax,1);
%loop over all the tdc_files 
for ii=1:iimax
    %predict the labview master trig time
    %use the write time to handle being unpacked from 7z
    tmp_est_labview_start=data.mcp_tdc.time_create_write(ii,2)...
        -anal_opts.trig_dld-anal_opts.dld_aquire-mean_delay_labview_tdc;
    [tval,nearest_idx]=closest_value(data.labview.time...
        ,tmp_est_labview_start);
    if abs(tval-tmp_est_labview_start)<time_thresh
        data.mcp_tdc.labview_shot_num(ii)=data.labview.shot_num(nearest_idx);
        data.mcp_tdc.probe.calibration(ii)=data.labview.calibration(nearest_idx);
    end 
end
clear('tmp_est_labview_start')

%% IMPORT THE ANALOG INPUT LOG
%load in the analog input files to check if the laser is single mode & that the potodiode value is close to the set
%point
% a two teired cache system is used one level for importing all 

anal_opts.ai_log.dir=anal_opts.tdc_import.dir;
anal_opts.ai_log.force_reimport=false;
anal_opts.ai_log.force_load_save=false;
anal_opts.ai_log.log_name='log_analog_in_';
%because im only passing the ai_log feild to aviod conflicts forcing a reimport i need to coppy these feilds
anal_opts.ai_log.calibration=data.mcp_tdc.probe.calibration;
anal_opts.ai_log.pd.set_probe=anal_opts.probe_set_pt;
anal_opts.ai_log.trig_dld=anal_opts.trig_dld;
anal_opts.ai_log.dld_aquire=anal_opts.dld_aquire;
anal_opts.ai_log.aquire_time=anal_opts.dld_aquire;
anal_opts.ai_log.trig_ai_in=anal_opts.trig_ai_in;
% set time matching conditions
anal_opts.ai_log.aquire_time=4;
anal_opts.ai_log.pd.diff_thresh=0.05;
anal_opts.ai_log.pd.std_thresh=0.05;
anal_opts.ai_log.pd.time_start=0.2;
anal_opts.ai_log.pd.time_stop=2;
anal_opts.ai_log.time_match_valid=8; %how close the predicted start of the shot is to the actual
%sfp options
anal_opts.ai_log.scan_time=1/20; %fast setting 1/100hz %estimate of the sfp scan time,used to set the window and the smoothing
anal_opts.ai_log.sfp.num_checks=inf; %how many places to check that the laser is single mode, inf=all scans
anal_opts.ai_log.sfp.peak_thresh=[-0.005,-0.005];%[0,-0.008]*1e-3; %theshold on the compressed signal to be considered a peak
anal_opts.ai_log.sfp.pzt_dist_sm=4.5;%minimum (min peak difference)between peaks for the laser to be considered single mode
anal_opts.ai_log.sfp.pzt_peak_width=0.15; %peak with in pzt voltage used to check that peaks are acually different and not just noise
anal_opts.ai_log.plot.all=false;
anal_opts.ai_log.plot.failed=true;

%do the ac waveform fit
anal_opts.ai_log.do_ac_mains_fit=false;

% Call the function
data.ai_log=ai_log_import(anal_opts.ai_log,data);

if isnan(anal_opts.probe_set_pt)
    stfig('finding pd setpt','add_stack',1);
    clf
    %plot the mean vs the pd std to determine what the setpt was
    plot(data.ai_log.pd.mean,data.ai_log.pd.std,'x')
    xlabel('mean pd voltage (v)')
    ylabel('std pd voltage(v)')
    %get some reasonable estimate for what the pd setpt was if it is unknown
    is_non_zero_mask=data.ai_log.pd.mean>0.1;
    %go one sd down from the mean pd variation during the probe
    std_upper_lim=mean(data.ai_log.pd.std(is_non_zero_mask))-std(data.ai_log.pd.std(is_non_zero_mask));
    std_upper_lim_mask=std_upper_lim<data.ai_log.pd.std;
    %then find the median value
    estimated_pd_setpt=median(data.ai_log.pd.mean(std_upper_lim_mask));
    yl=ylim;
    xl=xlim;
    hold on
    line([1,1]*estimated_pd_setpt,yl,'Color','k','LineWidth',1)
    %line([xl(1),xl(2)],[1,1]*(ai_log_single_out.pd.mean+ai_log_single_out.pd.std),'Color','r','LineWidth',3)
    %line([xl(1),xl(2)],[1,1]*(ai_log_single_out.pd.mean-ai_log_single_out.pd.std),'Color','r','LineWidth',3)
    
    anal_opts.ai_log.pd.set_probe=estimated_pd_setpt;
    data.ai_log=[];
    drawnow
    data.ai_log=ai_log_import(anal_opts.ai_log,data);
    
    %%write this out to the cal file
end

%%
%HACK IF SFP BROKEN
% 
% data.ai_log.ok.reg_pd=true(size(data.mcp_tdc.shot_num))';
% data.ai_log.ok.sfp=true(size(data.mcp_tdc.shot_num))';

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


%TODO: only pass anal_opts.wm_log

anal_opts.wm_log.plot.all=false;
anal_opts.wm_log.plot.failed=false;
anal_opts.wm_log.force_reimport=false;

anal_opts.wm_log.time_pd_padding=4; %check this many s each side of probe
anal_opts.wm_log.time_blue_padding=1; %check this many seconde each side of probe
anal_opts.wm_log.time_probe=3;
anal_opts.wm_log.ecd_volt_thresh=0.5;

anal_opts.wm_log.red_sd_thresh=5; %allowable standard deviation in MHz
anal_opts.wm_log.red_range_thresh=10; %allowable range deviation in MHz
anal_opts.wm_log.rvb_thresh=20; %allowable value of abs(2*red-blue)

anal_opts.wm_log.global=anal_opts.global;

data.wm_log.proc=wm_log_process(anal_opts,data);
clear('sub_data')

%TODO
% doublecheck that setpoint agrees with wavemeter value and isnt out by one shot
% clf
% plot(data.labview.setpoint-data.labview.setpoint(1),'xb')
% hold on
% plot(data.wm_log.proc.probe.freq.act.mean*1e6-data.labview.setpoint(1),'rx')




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
  

stfig('diagnostics and veto','add_stack',1);
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
legend('probe','all')
ylabel('Good?')
xlabel('Shot Number')
set(gca,'ytick',[0,1],'yticklabel',{'False','True'})
title('ALL ok')
ylim([-0.1,1.1])
tmp_num_shots=numel(data.mcp_tdc.shot_num);
tmp_num_ok_shots=sum(tmp_all_ok);
data.mcp_tdc.all_ok=tmp_all_ok;
data.mcp_tdc.probe.ok.all=tmp_probe_ok;
drawnow
fprintf('ok logic gives %u / %u shots for yeild %04.1f %%\n',...
    tmp_num_ok_shots,tmp_num_shots,1e2*tmp_num_ok_shots/tmp_num_shots)
% set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1600, 900])
% plot_name='check_logics';
% saveas(gcf,[anal_out.dir,plot_name,'.png'])
% saveas(gcf,[anal_out.dir,plot_name,'.fig'])


%% BINNING UP THE ATOM LASER PULSES
%now find the mean position of each pulse of the atom laser in each shot

%TODO: 
% some kind of guard pulses
% [x] a check that the pulse time seems reasonable
% convert the postions into velocity leaving the trap
anal_opts.atom_laser.pulsedt=8.000e-3;
anal_opts.atom_laser.start_pulse=1; %atom laser pulse to start with
anal_opts.atom_laser.pulses=136;
anal_opts.atom_laser.appr_osc_freq_guess=[52,40,40];
anal_opts.atom_laser.pulse_twindow=anal_opts.atom_laser.pulsedt*0.9;
anal_opts.atom_laser.xylim=anal_opts.tdc_import.txylim(2:3,:); %set same lims for pulses as import

anal_opts.atom_laser.plot.all=false;
%anal_opts.atom_laser.t0=0.417770;
anal_opts.atom_laser.global=anal_opts.global; %coppy global into the options structure
data.mcp_tdc.al_pulses=bin_al_pulses(anal_opts.atom_laser,data);


%% FITTING THE ATOM NUMBER
%use the inital few atom laser pulses in order to determine the atom number
%not of that much benifit TBH
anal_opts.atom_num_fit=[];
anal_opts.atom_num_fit.pulses=[1,20]; %min,max index of pulses
anal_opts.atom_num_fit.plot.each_shot=false;
anal_opts.atom_num_fit.plot.history=false;
anal_opts.atom_num_fit.qe=anal_opts.global.qe;

data.num_fit=fit_atom_number(anal_opts.atom_num_fit,data);


%% Load saved state
%save('before_fit.mat')
%load('before_fit.mat') % DEV DEV DEV

%% Correlate AC mains
% anal_opts.cancel_mains_corr.do=true;
% %anal_opts.cancel_mains_corr.tlim=[-3,2];
% %anal_opts.cancel_mains_corr.tsamp=1e-3;
% anal_opts.cancel_mains_corr.tlim=[-4,-1];
% anal_opts.cancel_mains_corr.tsamp=5e-4;
% data.corr_cancel=corr_ac_mains(data,anal_opts.cancel_mains_corr)

%% FITTING THE TRAP FREQUENCY

anal_opts.osc_fit.adaptive_freq=true; %estimate the starting trap freq 
anal_opts.osc_fit.dimesion=2; %Select coordinate to bin. 1=X, 2=Y.
anal_opts.osc_fit.appr_osc_freq_guess=[52,47.9,40];
anal_opts.osc_fit.freq_fit_tolerance=2; %hz arround the median to cut away
anal_opts.osc_fit.plot_fits=false;
anal_opts.osc_fit.plot_err_history=true;
anal_opts.osc_fit.plot_fit_corr=true;

anal_opts.osc_fit.global=anal_opts.global;
%%
data.osc_fit=fit_trap_freq(anal_opts.osc_fit,data);
%%
%data.osc_fit=fit_trap_freq_dev(anal_opts.osc_fit,data);

%% undo the aliasing
%this may need to change if the sampling freq changes
%initialize
data.osc_fit.trap_freq_recons=nan*data.osc_fit.ok.did_fits;
data.osc_fit.trap_freq_recons_unc=data.osc_fit.trap_freq_recons;
mask=data.osc_fit.ok.all; %set the masked values
data.osc_fit.trap_freq_recons(mask)=3*(1/anal_opts.atom_laser.pulsedt)+data.osc_fit.model_coefs(mask,2,1);
data.osc_fit.trap_freq_recons_unc(mask)=data.osc_fit.model_coefs(mask,2,2);

%% CHECK IF ATOM NUMBER DEPENDS ON PROBE BEAM
%figure
%clf
%plot(data.mcp_tdc.probe.freq.act.mean(data.osc_fit.ok.rmse),data.mcp_tdc.num_counts(data.osc_fit.ok.rmse),'x')

%% just plot with index
%plot(data.osc_fit.model_coefs(data.osc_fit.ok.rmse,2,1)'...
%    -mean(data.osc_fit.model_coefs(data.osc_fit.ok.rmse,2,1)))


%% plot not calibrations
%clf
%temp_cal=data.mcp_tdc.probe.calibration';
%temp_cal(isnan(temp_cal))=1;
%mask=data.osc_fit.ok.rmse & ~temp_cal & ~isnan(data.mcp_tdc.probe.freq.act.mean');
%plot(data.mcp_tdc.probe.freq.act.mean(mask),data.osc_fit.model_coefs(mask,2,1),'x')


%% create a model of the underlying trap frequency from the calibrations
% TODO
% - move signal calulation (trap freq difference) to this function
% - try other smoothing approaches
% - estimate error in difference/signal 
anal_opts.cal_mdl.smooth_time=60;
anal_opts.cal_mdl.plot=true;
anal_opts.cal_mdl.global=anal_opts.global;
data.cal=make_cal_model(anal_opts.cal_mdl,data);


%% calculae the probe beam trap frequency squared
% TODO
% - anharmonic correction
anal_opts.calc_sig=[];

data.signal=calculate_signal(anal_opts.calc_sig,data);


%% calculate blue probe freq
% convert freq to blue in hz apply aom shift to probe beam
data.blue_probe=calc_probe_blue(data.wm_log.proc,anal_opts.aom_freq);

%% Fit the tune out uisng all the data points from this scan
% fit using both a linear and quadratic polynomial
% TODO: optional fit weighting


anal_opts.fit_to_all=[];
anal_opts.fit_to_all.plot_inital=true;
anal_opts.fit_to_all.bootstrap=true;
%thresholds for CI
%sd         CI
%1          0.3174
%2          0.05
%3          2.699e-03
anal_opts.fit_to_all.sigma_disp=1;%one sd %confidence interval to display
anal_opts.fit_to_all.global=anal_opts.global;
anal_opts.fit_to_all.sigma_cut_outliers=3; %confidence interval for cutting outliers
anal_opts.fit_to_all.scale_x=1e-9;

data.to_fit_all=fit_to_all(anal_opts.fit_to_all,data);


%% segmented TO
%look at the tune out when fit to short segments
% TO DO, would be better if this called the fit_to script multiple times

anal_opts.fit_to_seg=[];
anal_opts.fit_to_seg.bootstrap=false;
anal_opts.fit_to_seg.plots=true;
anal_opts.fit_to_seg.clear_plot=true;
%thresholds for CI
%sd         CI
%1          1.3174
%2          0.05
%3          2.699e-03
anal_opts.fit_to_seg.sigma_disp=1;%one sd %confidence interval to display
anal_opts.fit_to_seg.global=anal_opts.global;
anal_opts.fit_to_seg.sigma_cut_outliers=3; %confidence interval for cutting outliers

anal_opts.fit_to_seg.scale_x=1e-9;
anal_opts.fit_to_seg.min_pts=7;

data.to_fit_seg=scan_segmented_fit_to(anal_opts.fit_to_seg,data);

%% write out the results

disp_to_results(data,anal_opts)

diary off 

%% what the output should have
% the aom offset that was included
% the directory that was used
% the results of the fit to all the data
% the results of the fit to each scan


%% damping results
%plot out what the distibution over damping times is
% figure(7);
% set(gcf,'color','w')
% histogram(1./data.osc_fit.model_coefs(data.osc_fit.ok.rmse,7,1),linspace(0,3,1e2))


%%
to_fit_seg=data.to_fit_seg;
to_fit_all=data.to_fit_all;

fprintf('saving results\n')
save(fullfile(anal_out.dir,'data_results.mat'),'to_fit_seg','to_fit_all','anal_opts','-nocompression','-v7.3')

fprintf('saving full output...')
%no compression bc its very slowwww
save(fullfile(anal_out.dir,'data_anal_full.mat'),'data','anal_opts','-nocompression','-v7.3')

%write a file called done to the out directory
fid = fopen(fullfile(anal_out.dir,'done.txt'),'wt');
fprintf(fid, 'done\n');
fclose(fid);
toc(main_trap_freq_timer)

fprintf('Done\n')

%catch err
%fprintf('Analysis on folder\n (%s) failed \n',anal_opts.tdc_import.dir) %Indicate if a directory couldn't be analysed properly
%msgText = getReport(err)
%end
catch e
    fprintf('caught error:%s',getReport(e))
    diary off
end %error catchign
end %process folder ?
end %loop over folders