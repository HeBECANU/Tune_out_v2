%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%interface between labview and matlab

%INPUT FROM LABVIEW
%   i          - integer,itteration number
%   auto_enable- boolean
%   file       - string : tells program where to look for setting files
%   file_exact - string
%   mloop      - boolean(legacy) do you want to interface with mloop or use matlab to scan over variables
%   control    - boolean(legacy)

%return
%exact_line- double (legacy just return a zero)
%new_path 

%to do
%convert to function
%Add conditional clauses to check for errors, write to log
%Add conditional to enable/disable settings updates - without having to restart LV? I guess one can
%just revert to v5 and reboot
%Move num_instr to instructions.m to put all customization in one file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%exact_line=0; %legacy to stop crashing

%setpoints = 362867621-6000; %debug only
%setpoints = linspace(362848466.40,362871075.28,50);
%setpoints = linspace(362853060,362880320,50); 
% setpoints = linspace(362868100-5000,362868100+5000,20); %to ~ 362868100(40)
%setpoints = linspace(362848466.40,362871075.28,50);
%setpoints = linspace(362853060,362880320,50); 
%freq_cen=(700939247.242651/2)-145.8;%145.9 %first big run of forbidden transtion %value from nist levels 
%freq_cen=(700939270.97/2)-151.95 ; %drake calc http://www.nrcresearchpress.com/doi/pdf/10.1139/p06-009
%freq_cen=350469637.449212;
freq_cen=362867621;
freq_range=2000;
freq_step=300;
setpoints = linspace(freq_cen-freq_range,freq_cen+freq_range,round(2*freq_range/freq_step)); 

path_log='Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\log_LabviewMatlab.txt';

%% Sequential setpoint
%avoid setting pointer to zero, i=iteration number passed from LabView
%add some random noise in the pointer to average out drifts a bit %round((rand(1)-0.5)*10
calibrate_interval=2;

if mod(i,calibrate_interval)==0
   %new_path='c:\remote\settings201826Nov224020.xml'; %probe off Lieutenant Angler
    %new_path='c:\remote\settings201826Nov224020.xml'; %probe off ,Blunderous Minstrel
    %new_path='c:\remote\settings201902Jan044412.xml';% set point 0.0v  '''Rich Carpenter'''
    %new_path='c:\remote\settings201904Jan004043.xml';% set point 0.0v  '''poor carpenter'''
    new_path='c:\remote\settings201917Jan112143.xml'; %probe off RF SW & flipper/shutt '''bourgeois  carpenter'''
    %write log entry
    f_log=fopen(path_log,'a');  % append to log-file
    nowdt=datetime('now');
    fprintf(f_log,'%.3f,%s,interfacev8,calibrate,itt,%u\n',posixtime(nowdt),datestr(nowdt,'yyyy-mm-ddTHH:MM:SS.FFF'),i);
    fclose(f_log);
else
    %new_path='c:\remote\settings201819Nov164232.xml'; %5v set point,Lieutenant Angler
    %new_path='c:\remote\settings201831Dec150742.xml';  %3v ,Blunderous Minstrel
    %new_path='c:\remote\settings201902Jan044348.xml';% set point 1.5v  '''Rich Carpenter'''
    %new_path='c:\remote\settings201904Jan004030.xml'% set point 0.0v  '''poor carpenter'''
    %new_path='c:\remote\settings201915Jan050928.xml';% set point 1.25V '''bourgeois carpenter'''
    %new_path='c:\remote\settings201901Feb020412.xml';% set point 0.31V '''bourgeois carpenter'''
    %new_path='c:\remote\settings201909Feb010057.xml';% set point 1.00V '''bourgeois carpenter'''
    %new_path='c:\remote\settings201912Feb171618.xml';% set point 0.4V '''bourgeois carpenter'''
    %new_path='c:\remote\settings201913Feb095900.xml';% set point 0.8V '''bourgeois carpenter'''
    %new_path='c:\remote\settings201913Feb124058.xml';% set point 0.7V '''bourgeois carpenter'''
    %new_path='c:\remote\settings201913Feb154835.xml';% set point 0.9V '''bourgeois carpenter'''
    %new_path='c:\remote\settings201913Feb201434.xml';% set point 0.4V '''bourgeois carpenter'''
    %new_path='c:\remote\settings201914Feb121756.xml'% set point 0.6V '''bourgeois carpenter'''
    %new_path='c:\remote\settings201914Feb162456.xml';% set point 2.0V '''bourgeois carpenter'''
    %new_path='c:\remote\settings201914Feb181533.xml';% set point 2.5V '''bourgeois carpenter'''
    new_path='c:\remote\settings201914Feb193925.xml';% set point 1.0V '''bourgeois carpenter'''
    
    %new_path='c:\remote\settings201909Feb202659.xml';% set point 0.20V '''bourgeois carpenter'''
    pointer=floor(i/calibrate_interval)*(calibrate_interval-1)+ rem(i,calibrate_interval);
    pointer = mod(pointer-1,length(setpoints))+1; 
    setpt = setpoints(pointer);
%     
    %if there is no wm feedback running this blocks execution
%      t = tcpip('0.0.0.0', 33333, 'NetworkRole', 'server');
%      fopen(t)
%      fwrite(t, setpt,'double')
%      fclose(t)
%      
    %write log entry
    f_log=fopen(path_log,'a');  % append to log-file
    nowdt=datetime('now');
    fprintf(f_log,'%.3f,%s,interfacev8,measure_probe,%f,itt,%u\n',posixtime(nowdt),datestr(nowdt,'yyyy-mm-ddTHH:MM:SS.FFF'),setpt,i);
    fclose(f_log);
end


%code scraps
%dir_this=fileparts(mfilename('fullpath'));   % get dir of this script
%path_log=strcat(dir_this,'\logs\','log_LabviewMatlab.txt');    % path of log file




