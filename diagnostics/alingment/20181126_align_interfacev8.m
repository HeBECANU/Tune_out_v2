
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%DOCUMENTATION
%bryce:what needs to be in the path for this to run properly
%explain briefly what this program does and what the user will need to edit

%to do
%could read in number of itterations from settings file in order to be a
%bit inelegent
%Add conditional clauses to check for errors, write to log
%Add conditional to enable/disable settings updates - without having to restart LV? I guess one can
%just revert to v5 and reboot
%Move num_instr to instructions.m to put all customization in one file

%INPUT FROM LABVIEW
%mloop,boolean : do you want to interface with mloop or use matlab to scan over variables
%file,string : tells program where to look for setting files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%setpoints = linspace(362848466.40,362871075.28,50);
%setpoints = linspace(362853060,362880320,50); 
setpoints = linspace(362868100-6000,362868100+6000,20); %to ~ 362868100(40)

    
path_user_config='..\MloopUserConfig.txt';  % this is the user config file regardless of Mloop
dir_this=fileparts(mfilename('fullpath'));   % get dir of this script
%path_log=strcat(dir_this,'\logs\','log_LabviewMatlab.txt');    % path of log file
path_log='Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\log_LabviewMatlab.txt';

%% Sequential setpoint
%avoid setting pointer to zero, i=iteration number passed from LabView
%add some random noise in the pointer to average out drifts a bit %round((rand(1)-0.5)*10
calibrate_interval=2;


if mod(i,calibrate_interval)==0
    %new_path='c:\remote\settings201819Nov193558.xml'; %probe off Lieutenant Angler
    new_path='c:\remote\settings201826Nov224020.xml';
    %write log entry
    f_log=fopen(path_log,'a');  % append to log-file
    nowdt=datetime('now');
    fprintf(f_log,'%.3f,%s,interfacev8,calibrate,itt,%u\n',posixtime(nowdt),datestr(nowdt,'yyyy-mm-ddTHH:MM:SS.FFF'),i)
    fclose(f_log);
else
    %new_path='c:\remote\settings201819Nov164232.xml'; %5v set point,Lieutenant Angler
    new_path='c:\remote\settings201826Nov212303.xml';  %5v ,Lieutenant  daggertooth
    pointer=floor(i/calibrate_interval)*(calibrate_interval-1)+ rem(i,calibrate_interval);
    pointer = mod(pointer-1,length(setpoints))+1; 
    setpt = setpoints(pointer);
%     
%     %if there is no wm feedback running this blocks execution
%      t = tcpip('0.0.0.0', 33333, 'NetworkRole', 'server');
%      fopen(t)
%      fwrite(t, setpt,'double')
%      fclose(t)
    
    %write log entry
    f_log=fopen(path_log,'a');  % append to log-file
    nowdt=datetime('now');
    fprintf(f_log,'%.3f,%s,interfacev8,measure_probe,%f,itt,%u\n',posixtime(nowdt),datestr(nowdt,'yyyy-mm-ddTHH:MM:SS.FFF'),setpt,i)
    fclose(f_log);
end



%% Random setpoint
%setpt = datasample(setpoints,1); %Picks a setpoint at random

%% Scanning set point
% num_shots = 10;
% pointer = mod(idivide(int32(i),int32(num_shots),length(setpoints));
% setpt = setpoints(pointer);


% t = tcpip('0.0.0.0', 33333, 'NetworkRole', 'server');
% fopen(t)
% fwrite(t, 362871075.28,'double')
% fclose(t)



%new_path = file; %Retain labview settings




