function back_calc_out=background_single_file(args_single)

%------------- BEGIN USER VAR --------------
args_single.pzt_attenuation_factor=0.2498; %set by the voltage divider box between pzt driver and daq
%------------- END USER VAR --------------

%------------- BEGIN CODE --------------






if args_single.dir(end) ~= filesep, args_single.dir = [args_single.dir, filesep]; end


%%load the data
path=fullfile(args_single.dir,args_single.fname);
fid = fopen(path,'r');
raw_line = fread(fid,Inf,'*char')'; % this is the fastest string read method 3.5s for 20 files
fclose(fid);
ai_dat=jsondecode(raw_line);
samples= size(ai_dat.Data,2);
sr=ai_dat.sample_rate;
aquire_time=samples/sr;

ai_dat.time=(0:(samples-1))/sr;%not sure if this is off by 1 sample in time


if args_single.plot.all    
    stfig('analog_ins','add_stack',1);
    clf;
    plot(ai_dat.time,ai_dat.Data,'b')
    ylabel('voltage')
    xlabel('time (s)')
    pause(1e-6)

end

%% Process the scan to find the background power as a function of integration width
% dev load('ai_log_state_before_is_sm.mat')
%load('ai_log_state_before_is_sm.mat') %DEV DEV DEV REMOVE FOR USE

%%TODO
% input label check
sfp_pzt= col_vec(ai_dat.Data(3,:));
sfp_pd_raw= col_vec(ai_dat.Data(2,:));
sfp_pd_cmp= col_vec(ai_dat.Data(4,:));

%correct for the attenuation of the pzt voltage that goes into the DAQ
sfp_pzt=sfp_pzt/args_single.pzt_attenuation_factor;


% set up input struct for is_laser_single_mode
back_calc_in.pd_voltage=[sfp_pd_raw,sfp_pd_cmp];
back_calc_in.times= col_vec(ai_dat.time);
back_calc_in.pzt_voltage= sfp_pzt;

%move the below options back to main trap freq
back_calc_in.scan_type='sawtooth';
back_calc_in.pzt_scan_period=1;  %estimate of the sfp scan time,used to set the window and the smoothing
back_calc_in.pd_filt_factor=1e-3; %fraction of a scan to smooth the pd data by for peak detection
back_calc_in.ptz_filt_factor_pks=1e-3;  %fraction of a scan to smooth the pzt data by for peak detection
back_calc_in.pzt_filt_factor_deriv=1e-3; %fraction of a scan to smooth the data by for derivative detection
back_calc_in.sfp_finesse=args_single.sfp_finesse;
back_calc_in.int_width=args_single.int_width;
back_calc_in.pd_amp_min=1; %minimum range of the pd signal to indicate the laser has sufficient power
back_calc_in.skip_wf_check=0;
back_calc_in.plot=args_single.plot;

back_calc_out=background_calculate(back_calc_in);


end