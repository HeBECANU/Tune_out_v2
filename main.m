%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine the Tune out from a dataset
% using the measured change in trap frequency from the application of a probe beam.
% application of the tune out probe beam.
% The script:
%	* defines the user controled options
%   * Imports all the tdc data files 
%	* Imports labview log file
%   * Imports the wavemeter log file
%   * Match upt the Labview data withthe tdc data
%   * Imports the analog in log file , for each file the import:
%       * checks that the pd voltage is ok
%       * check that the laser is single mode using the scanning fabry perot signals
%   * Check that the wavemeter readings are ok for each shot
%       *checks that the wavelengths is stable during the probe intterogation
%       *checks that the red wavelength is ~half the blue
%       *checks that the double pd voltage is ok (now redundant beacuse of probe pd)
%   * checks that the number of counts in the file is ok
%   * combines all these checks into one master check
%   * bins up each pulse of the AL
%   * Fits the trap frequency
%	* Investigate fit correlations
%	* mask out only the (good)calibrations shots and make a model of how the (unpeturbed) trap freq changes in time
%   * plot out (non calibration) (good) data and then fit the probe beam
%       wavelength, identifying the tuneout wavelength and giving a
%       stastistical uncertainty.

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
%   - maybe process the analog logs into a output file for faster reading in??
%   - weight the final TO fit by the trap freq fit unc
%	-make unique plot numbers
%   -the fit error depends on wavelength indicating that the model does not
%   have enough freedom
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
% Last revision:2018-10-09


%close all
clear all
tic
%%
% BEGIN USER VAR-------------------------------------------------
anal_opts=[];
anal_opts.tdc_import.dir='Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20190115_baseline_to_1\';
anal_opts.tdc_import.file_name='d';
anal_opts.tdc_import.force_load_save=false;   %takes precidence over force_reimport
anal_opts.tdc_import.force_reimport=false;
anal_opts.tdc_import.force_forc=false;
anal_opts.tdc_import.dld_xy_rot=0.61;
%Should probably try optimizing these
tmp_xlim=[-30e-3, 30e-3];     %tight XY lims to eliminate hot spot from destroying pulse widths
tmp_ylim=[-30e-3, 30e-3];
tlim=[0,4];
anal_opts.tdc_import.txylim=[tlim;tmp_xlim;tmp_ylim];

anal_opts.max_runtime=inf;%cut off the data run after some number of hours
anal_opts.atom_laser.pulsedt=8.000e-3;
anal_opts.atom_laser.t0=0.41784; %center i ntime of the first pulse
anal_opts.atom_laser.start_pulse=1; %atom laser pulse to start with
anal_opts.atom_laser.pulses=100;
anal_opts.atom_laser.appr_osc_freq_guess=[52,40,40];
anal_opts.atom_laser.pulse_twindow=anal_opts.atom_laser.pulsedt*0.9;
anal_opts.atom_laser.xylim=anal_opts.tdc_import.txylim(2:3,:); %set same lims for pulses as import

anal_opts.global.fall_time=0.417;
anal_opts.global.qe=0.09;

anal_opts.trig_dld=20.3;
anal_opts.dld_aquire=4;
anal_opts.aquire_time=4;
anal_opts.trig_ai_in=20;


anal_opts.wm_log.plot_all=true;
anal_opts.wm_log.plot_failed=true;


anal_opts.osc_fit.binsx=1000;
anal_opts.osc_fit.blur=1;
anal_opts.osc_fit.xlim=[-20,20]*1e-3;
anal_opts.osc_fit.tlim=[0.86,1.08];
anal_opts.osc_fit.dimesion=2; %Select coordinate to bin. 1=X, 2=Y.

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
%% IMPORT TDC DATA to data.mcp_tdc
anal_opts.tdc_import.shot_num=find_data_files(anal_opts.tdc_import);
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



%% Match Labview data
%because the mcp-dld detector has no direct communication with the bec computer
% the data.labview.shot_num does not nessesarily correspond to data.mcp_tdc.shot_num


    %try and match up the file with if it is a calibaration using the time
    %it is slightly overkill here to search each one, but just being extra
    %cautious/flexible
    time_thresh=4; %how close for the times to be considered the same shot
    %lets examine what the time difference does
    sfigure(45);
    set(gcf,'color','w')
    clf
    imax=min([size(data.labview.time,2),size(data.mcp_tdc.time_create_write,1)]);
    %imax=5000;
    time_diff=data.mcp_tdc.time_create_write(1:imax,2)'-anal_opts.dld_aquire-anal_opts.trig_dld-...
        data.labview.time(1:imax);
    mean_delay_labview_tdc=mean(time_diff);
    plot(time_diff)
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
            -anal_opts.trig_dld-anal_opts.dld_aquire-mean_delay_labview_tdc;
        [tval,nearest_idx]=closest_value(data.labview.time...
            ,est_labview_start);
        nearest_idx
        if abs(tval-est_labview_start)<time_thresh
            data.mcp_tdc.labview_shot_num(ii)=data.labview.shot_num(nearest_idx);
            data.mcp_tdc.probe.calibration(ii)=data.labview.calibration(nearest_idx);
        end 
    end

%% IMPORT THE ANALOG INPUT LOG
%the code will check that the probe beam PD was ok and that the laser was single mode
anal_opts.ai_log.dir=anal_opts.tdc_import.dir;
anal_opts.ai_log.force_reimport=false;
anal_opts.ai_log.force_load_save=false;
anal_opts.ai_log.log_name='log_analog_in_';
anal_opts.ai_log.pd.set=data.mcp_tdc.probe.calibration;
%nan compatable logical inverse
anal_opts.ai_log.pd.set(~isnan(anal_opts.ai_log.pd.set))=~anal_opts.ai_log.pd.set(~isnan(anal_opts.ai_log.pd.set))
anal_opts.ai_log.pd.set=anal_opts.ai_log.pd.set*3;
%anal_opts.ai_log.pd.set(isnan(anal_opts.ai_log.pd.set))=0;
anal_opts.ai_log.aquire_time=4;
anal_opts.ai_log.pd.diff_thresh=0.1;
anal_opts.ai_log.pd.std_thresh=0.1;
anal_opts.ai_log.pd.time_start=0.2;
anal_opts.ai_log.pd.time_stop=2;
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

%% Call the function
ai_log_out=ai_log_import(anal_opts.ai_log,data);
%copy the output across
data.ai_log=ai_log_out;

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
anal_opts.wm_log.plot_failed=true;
anal_opts.wm_log.force_reimport=false;

anal_opts.wm_log.time_pd_padding=4; %check this many s each side of probe
anal_opts.wm_log.time_blue_padding=1; %check this many seconde each side of probe
anal_opts.wm_log.time_probe=3;
anal_opts.wm_log.ecd_volt_thresh=0.5;

anal_opts.wm_log.red_sd_thresh=50; %allowable standard deviation in MHz
anal_opts.wm_log.red_range_thresh=50; %allowable range deviation in MHz
anal_opts.wm_log.rvb_thresh=10; %allowable value of abs(2*red-blue)

data.wm_log.proc=wm_log_process(anal_opts,data);
clear('sub_data')


%% CHECK ATOM NUMBER
sfigure(1)
clf
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
  

figure(12);
clf
set(gcf,'color','w')
subplot(2,1,1)
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
subplot(2,1,2)
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
data.mcp_tdc.all_ok=tmp_all_ok;
data.mcp_tdc.probe.ok.all=tmp_probe_ok;


fprintf('ok logic gives %u / %u shots for yeild %04.1f %%\n',...
    tmp_num_ok_shots,tmp_num_shots,1e2*tmp_num_ok_shots/tmp_num_shots)
% set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1600, 900])
% plot_name='check_logics';
% saveas(gcf,[anal_out.dir,plot_name,'.png'])
% saveas(gcf,[anal_out.dir,plot_name,'.fig'])

%% BINNING UP THE ATOM LASER PULSES
%now find the mean position of each pulse of the atom laser in each shot
data.mcp_tdc.al_pulses=bin_al_pulses(anal_opts.atom_laser,data);

%% FITTING THE TRAP FREQUENCY
anal_opts.osc_fit.adaptive_freq=true; %estimate the starting trap freq 
anal_opts.osc_fit.appr_osc_freq_guess=[52,47.9,40];
anal_opts.osc_fit.freq_fit_tolerance=2; %hz arround the median to cut away
anal_opts.osc_fit.plot_fits=true;
anal_opts.osc_fit.plot_err_history=true;
anal_opts.osc_fit.plot_fit_corr=true;

anal_opts.osc_fit.global=anal_opts.global;
data.osc_fit=fit_trap_freq(anal_opts.osc_fit,data);

%% undo the aliasing
%this may need to change if the sampling freq changes

data.osc_fit.trap_freq_recons=nan*data.osc_fit.ok.did_fits;
mask=data.osc_fit.ok.all;
data.osc_fit.trap_freq_recons(mask)=3*(1/anal_opts.atom_laser.pulsedt)+data.osc_fit.model_coefs(mask,2,1);


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
anal_opts.cal_mdl.smooth_time=100;
anal_opts.cal_mdl.plot=true;
anal_opts.cal_mdl.global=anal_opts.global;
data.cal=make_cal_model(anal_opts.cal_mdl,data);


%% segmented TO
%look at the tune out when fit to short segments
% TO DO, would be better if this called the fit_to script multiple times
anal_opts.fit_to=[];
anal_opts.fit_to.bootstrap=false;
anal_opts.fit_to.plots=false;
%thresholds for CI
%sd         CI
%1          0.3174
%2          0.05
%3          2.699e-03
anal_opts.fit_to.ci_size_disp=0.3174;%one sd %confidence interval to display
anal_opts.fit_to.global=anal_opts.global;
anal_opts.fit_to.ci_size_cut_outliers=0.05; %confidence interval for cutting outliers
anal_opts.fit_to.scale_x=1e-9;
anal_opts.fit_to.min_pts=10;

anal_opts.fit_to.seg_time=60*30;
anal_opts.fit_to.seg_shift=1*anal_opts.fit_to.seg_time;
to_seg_fits=segmentd_fit_to(anal_opts.fit_to,data);



%% Fit the Tune Out
% anal_opts.fit_to=[];
% anal_opts.fit_to.plot_inital=true;
% anal_opts.fit_to.bootstrap=true;
% %thresholds for CI
% %sd         CI
% %1          0.3174
% %2          0.05
% %3          2.699e-03
% anal_opts.fit_to.ci_size_disp=0.3174;%one sd %confidence interval to display
% anal_opts.fit_to.global=anal_opts.global;
% anal_opts.fit_to.ci_size_cut_outliers=0.05; %confidence interval for cutting outliers
% anal_opts.fit_to.scale_x=1e-9;
% 
% to_res=fit_to(anal_opts.fit_to,data);
% data.to_fit=to_res;
% 
% to_fit_trimed_val=to_res.fit_trimmed.to_freq;
% to_fit_unc_boot=to_res.fit_trimmed.to_unc_boot;
% to_fit_unc_fit=to_res.fit_trimmed.to_unc_fit;
% to_fit_unc_unc_boot=to_res.fit_trimmed.boot.se_se_opp/anal_opts.fit_to.scale_x;

%% Analyse the effect of nonlinear terms on Tune out
anal_opts.fit_to=[];
anal_opts.fit_to.plot_inital=true;
anal_opts.fit_to.bootstrap=true;
%thresholds for CI
%sd         CI
%1          0.3174
%2          0.05
%3          2.699e-03
anal_opts.fit_to.ci_size_disp=0.3174;%one sd %confidence interval to display
anal_opts.fit_to.global=anal_opts.global;
anal_opts.fit_to.ci_size_cut_outliers=0.05; %confidence interval for cutting outliers
anal_opts.fit_to.scale_x=1e-9;

to_res=fit_to_nonlin(anal_opts.fit_to,data);
data.to_fit_nonlin=to_res;

to_fit_trimed_val=to_res.fit_trimmed.to_freq;
to_fit_unc_boot=to_res.fit_trimmed.to_unc_boot;
to_fit_unc_fit=to_res.fit_trimmed.to_unc_fit;
to_fit_unc_unc_boot_lin=to_res.fit_trimmed.boot{1}.se_se_opp/anal_opts.fit_to.scale_x;
to_fit_unc_unc_boot_quad=to_res.fit_trimmed.boot{2}.se_se_opp/anal_opts.fit_to.scale_x;

%% write out the results
%inverse scaled gradient to give the single shot uncert (with scaling factor to include calibration)
tot_num_shots=to_res.num_shots+data.cal.num_shots;
single_shot_uncert=to_res.fit_trimmed.single_shot_uncert_boot{1}...
    *sqrt(tot_num_shots/to_res.num_shots);
fprintf('\n====TO fit results==========\n')
fprintf('dir =%s\n',anal_opts.tdc_import.dir)
fprintf('median damping time %.2f\n',median(1./data.osc_fit.model_coefs(data.osc_fit.ok.rmse,7,1)))
%calculate some statistics and convert the model parameter into zero crossing and error therin
old_to_wav=413.0938e-9;
new_to_freq_unc=to_fit_unc_boot;
%to_res.fit_trimmed.to_unc_fit
to_wav_val_lin=const.c/(to_fit_trimed_val{1}*2);
to_wav_unc_lin=new_to_freq_unc{1}*const.c/((to_fit_trimed_val{1}*2)^2);
to_wav_val_quad=const.c/(to_fit_trimed_val{2}*2);
to_wav_unc_quad=new_to_freq_unc{2}*const.c/((to_fit_trimed_val{2}*2)^2);
fprintf('run start time               %.1f (posix)\n',...
    data.mcp_tdc.time_create_write(1,2)-anal_opts.trig_dld-anal_opts.dld_aquire)
fprintf('run stop time                %.1f (posix)\n',...
    data.mcp_tdc.time_create_write(end,2)-anal_opts.trig_dld-anal_opts.dld_aquire)
fprintf('duration                     %.1f (s)\n',...
    data.mcp_tdc.time_create_write(end,2)-data.mcp_tdc.time_create_write(1,2))
fprintf('TO freq (Linear)             %.1f±(%.0f±%.0f) MHz\n',...
    to_fit_trimed_val{1}*1e-6,new_to_freq_unc{1}*1e-6,to_fit_unc_unc_boot_lin*1e-6)
fprintf('TO wavelength (Linear)       %.6f±%f nm \n',to_wav_val_lin*1e9,to_wav_unc_lin*1e9)
fprintf('TO freq (Quadratic)          %.1f±(%.0f±%.0f) MHz\n',...
    to_fit_trimed_val{2}*1e-6,new_to_freq_unc{2}*1e-6,to_fit_unc_unc_boot_quad*1e-6)
fprintf('TO wavelength (Quadratic)    %.6f±%f nm \n',to_wav_val_quad*1e9,to_wav_unc_quad*1e9)
fprintf('diff between Lin and Quad    %e±%e nm \n',(to_wav_val_lin-to_wav_val_quad)*1e9,sqrt(to_wav_unc_lin^2+to_wav_unc_quad^2)*1e9)
fprintf('diff from TOV1               %e±%e nm \n',(to_wav_val_lin-old_to_wav)*1e9,to_wav_unc_lin*1e9)
%more logic needs to be included here
fprintf('number of probe files        %u \n',to_res.num_shots)
fprintf('number of calibration files  %u \n',data.cal.num_shots)
fprintf('total used                   %u \n',tot_num_shots)
fprintf('files with enough number     %u\n',sum(data.mcp_tdc.num_ok'))
fprintf('shot uncert scaling @1SD %.1f MHz, %.2f fm /sqrt(shots)\n',single_shot_uncert*1e-6,...
    single_shot_uncert*const.c/((to_fit_trimed_val{1}*2)^2)*10^15)
%predicted uncert using this /sqrt(n), unless derived differently this is pointless
%fprintf('predicted stat. uncert %.1f MHz, %.2f fm\n',single_shot_uncert/sqrt(tot_num_shots)*1e-6,...
%    single_shot_uncert/sqrt(tot_num_shots)*const.c/((to_fit_trimed_val*2)^2)*10^15)

diary off


%%
fprintf('saving output...')
%no compression bc its very slowwww
%save(fullfile(anal_out.dir,'data_anal_full.mat'),'data','to_res','anal_opts','-nocompression','-v7.3')
fprintf('Done')

%% try and find what the outliers were doing
%given this mask ~color_idx find the shot nums and times 
%~isnan(probe_freq)


%% damping results
%plot out what the distibution over damping times is
% figure(7);
% set(gcf,'color','w')
% histogram(1./data.osc_fit.model_coefs(data.osc_fit.ok.rmse,7,1),linspace(0,3,1e2))

toc