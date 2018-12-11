% to determine the alingment we have attempted to use the 2^3 S_1 -> 2^3P_0 transtion with the curcular polarization
% that is not not absorbed when the atoms are perfectly aligned with the beam
% to determine the transtion rate i did a pulsed atom laser (with constant power) then the pulses were paused and the
% probe(1803nm) beam was turned on for some time afterwards the pulse atom laser was restarted
% will use two atom number fits to these PAL segments to determine the transtion rate
% critialy if it changes as the nuller bias is changed

clear all
tic
%%
% BEGIN USER VAR-------------------------------------------------
anal_opts=[];
%anal_opts.tdc_import.dir='Y:\EXPERIMENT-DATA\Tune Out V2\20180826_testing_wm_log\';
anal_opts.tdc_import.dir='Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\setting_up_1083_mag_align\20181210_nuller_-7.0_qt_nan_1083_off';
anal_opts.tdc_import.file_name='d';
anal_opts.tdc_import.force_load_save=false;   %takes precidence over force_reimport
anal_opts.tdc_import.force_reimport=false;
anal_opts.tdc_import.force_forc=false;
anal_opts.tdc_import.dld_xy_rot=0.61;

tmp_xlim=[-30e-3, 30e-3];     %tight XY lims to eliminate hot spot from destroying pulse widths
tmp_ylim=[-30e-3, 30e-3];
tlim=[0,4];
anal_opts.tdc_import.txylim=[tlim;tmp_xlim;tmp_ylim];

anal_opts.max_runtime=10;%cut off the data run after some number of hours
anal_opts.atom_laser.pulsedt=8.000e-3;
anal_opts.atom_laser.t0=0.41784; %center i ntime of the first pulse
anal_opts.atom_laser.start_pulse=1; %atom laser pulse to start with
anal_opts.atom_laser.pulses=150;
anal_opts.atom_laser.appr_osc_freq_guess=[52,40,40];
anal_opts.atom_laser.pulse_twindow=anal_opts.atom_laser.pulsedt*0.9;
anal_opts.atom_laser.xylim=anal_opts.tdc_import.txylim(2:3,:); %set same lims for pulses as import

anal_opts.global.fall_time=0.417;
anal_opts.global.qe=0.09;

anal_opts.trig_dld=20.3;
anal_opts.dld_aquire=4;
anal_opts.aquire_time=4;
anal_opts.trig_ai_in=20;
anal_opts.aom_freq=0;%190*1e6;%Hz %set to zero for comparison with previous data runs
anal_opts.probe_set_pt=3.0;

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

%% OK
data.mcp_tdc.all_ok=data.mcp_tdc.num_counts>10;

%% BINNING UP THE ATOM LASER PULSES
%now find the mean position of each pulse of the atom laser in each shot
data.mcp_tdc.al_pulses=bin_al_pulses(anal_opts.atom_laser,data);


%% FITTING THE ATOM NUMBER PART 1
%use the inital few atom laser pulses in order to determine the atom number
%not of that much benifit TBH
anal_opts.atom_num_fit=[];
anal_opts.atom_num_fit.pulses=[1,15]; %min,max index of pulses
anal_opts.atom_num_fit.plot.each_shot=false;
anal_opts.atom_num_fit.plot.history=false;
anal_opts.atom_num_fit.qe=anal_opts.global.qe;

data.num_fit.pre_probe=fit_atom_number(anal_opts.atom_num_fit,data);

%% FITTING THE ATOM NUMBER PART 2
anal_opts.atom_num_fit.pulses=[52,120]; %min,max index of pulses
anal_opts.atom_num_fit.plot.each_shot=false;
data.num_fit.post_probe=fit_atom_number(anal_opts.atom_num_fit,data);

%% Compare
 
atom_num_pre=cellfun(@(x) x(16),data.num_fit.pre_probe.fit_predict);
atom_num_post=cellfun(@(x) x(52),data.num_fit.post_probe.fit_predict);
frac_rem_vec_lin=atom_num_post./atom_num_pre;
frac_rem_vec_log=log(atom_num_post./atom_num_pre);
attn_log_mean=mean(frac_rem_vec_log);
attn_log_std=std(frac_rem_vec_log)/numel(frac_rem_vec_log);
attn_lin_mean=mean(frac_rem_vec_lin);
attn_lin_std=std(frac_rem_vec_lin)/numel(frac_rem_vec_lin);
fprintf('frac atoms after/before probe %.4f±%.4f \n',attn_lin_mean,attn_lin_std)
fprintf('log frac atoms after/before probe %.4f±%.4f \n',attn_log_mean,attn_log_std)


