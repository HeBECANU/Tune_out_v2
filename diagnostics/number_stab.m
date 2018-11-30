



anal_opts=[]
%anal_opts.tdc_import.dir='Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20181017_probe_off_trap_freq_drift\';
anal_opts.tdc_import.dir='Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\';
anal_opts.tdc_import.file_name='d';
anal_opts.tdc_import.force_load_save=false;   %takes precidence over force_reimport
anal_opts.tdc_import.force_reimport=false;
anal_opts.tdc_import.force_forc=false;
anal_opts.tdc_import.dld_xy_rot=0.61;

tmp_xlim=[-30e-3, 30e-3];     %tight XY lims to eliminate hot spot from destroying pulse widths
tmp_ylim=[-30e-3, 30e-3];
tlim=[0,4];
anal_opts.tdc_import.txylim=[tlim;tmp_xlim;tmp_ylim];




%go up a directory and then add all subfolders
%add all subfolders
folder = fileparts(which(mfilename));
folder=strsplit(folder,filesep); %go up a directory
folder=strjoin(folder(1:end-1),filesep);
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));



%%

hebec_constants
anal_opts.tdc_import.mat_save=false;

if anal_opts.tdc_import.dir(end) ~= '\', dirpath = [dirpath '\']; end
if (exist([anal_opts.tdc_import.dir,'out'], 'dir') == 0), mkdir([anal_opts.tdc_import.dir,'out']); end
 
anal_out.dir=sprintf('%sout\\monitor\\',...
    anal_opts.tdc_import.dir);
if (exist(anal_out.dir, 'dir') == 0), mkdir(anal_out.dir); end
anal_opts.global.out_dir=anal_out.dir;




%%

anal_opts.tdc_import.shot_num=find_data_files(anal_opts.tdc_import);
anal_opts.tdc_import.shot_num=anal_opts.tdc_import.shot_num(200:end);
%anal_opts.tdc_import.shot_num=47:65;

%max_shot_num=max(anal_opts.tdc_import.shot_num);
%anal_opts.tdc_import.shot_num=anal_opts.tdc_import.shot_num(...
%    anal_opts.tdc_import.shot_num>(max_shot_num-anal_opts.history.shots));

data.mcp_tdc=import_mcp_tdc_data(anal_opts.tdc_import);


%%
sfigure(5)
set(gcf,'color','w')
time_since_start=data.mcp_tdc.time_create_write(:,2)-data.mcp_tdc.time_create_write(1,2);
%JR set second index to 2 because of anomalous delay in create time (fixed in write time)
subplot(1,2,1)
plot(time_since_start,data.mcp_tdc.num_counts,'k')

% % JR playing with smoothed trace
% low_shots = data.mcp_tdc.num_counts < 200;
% gapped_mean = [0,0.5*(data.mcp_tdc.num_counts(3:end) + data.mcp_tdc.num_counts(1:end-2)),0];
% filled_x = data.mcp_tdc.num_counts.*(1-low_shots) + gapped_mean.*low_shots;
% plot(time_since_start,filled_x,'k')

xlabel('time(s)')
ylabel('detected counts')
subplot(1,2,2)
fft_out=fft_tx(time_since_start,data.mcp_tdc.num_counts,10);
plot(1./fft_out(1,:),abs(fft_out(2,:)),'k'); % plot as a function of the period
%plot(fft_out(1,:),abs(fft_out(2,:))); % plot as a function of the frequency
xlabel('period (s)')
ylabel('count amplitude')

