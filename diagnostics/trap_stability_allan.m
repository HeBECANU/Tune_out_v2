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
anal_opts=[]
anal_opts.tdc_import.dir='Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20181010_probe_off_trap_freq_drift\';
anal_opts.tdc_import.file_name='d';
anal_opts.tdc_import.force_load_save=false;   %takes precidence over force_reimport
anal_opts.tdc_import.force_reimport=false;
anal_opts.tdc_import.force_forc=false;
anal_opts.tdc_import.dld_xy_rot=0.61;

tmp_xlim=[-30e-3, 30e-3];     %tight XY lims to eliminate hot spot from destroying pulse widths
tmp_ylim=[-30e-3, 30e-3];
tlim=[0,4];
anal_opts.tdc_import.txylim=[tlim;tmp_xlim;tmp_ylim];


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
anal_opts.trig_ai_in=20;


anal_opts.osc_fit.binsx=1000;
anal_opts.osc_fit.blur=1;
anal_opts.osc_fit.xlim=[-20,20]*1e-3;
anal_opts.osc_fit.tlim=[0.86,1.08];
anal_opts.osc_fit.dimesion=2; %Select coordinate to bin. 1=X, 2=Y.

anal_opts.history.shots=200;

% END USER VAR-----------------------------------------------------------

%add all subfolders
folder = fileparts(which(mfilename));
folder=strsplit(folder,filesep); %go up a directory
folder=strjoin(folder(1:end-1),filesep);
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));

hebec_constants
anal_opts.tdc_import.mat_save=false;
anal_opts.global.velocity=const.g0*anal_opts.global.fall_time;

if anal_opts.tdc_import.dir(end) ~= '\', dirpath = [dirpath '\']; end
if (exist([anal_opts.tdc_import.dir,'out'], 'dir') == 0), mkdir([anal_opts.tdc_import.dir,'out']); end
 
anal_out.dir=sprintf('%sout\\monitor\\',...
    anal_opts.tdc_import.dir);
if (exist(anal_out.dir, 'dir') == 0), mkdir(anal_out.dir); end
anal_opts.global.out_dir=anal_out.dir;




%%
loop_num=0;
sfigure(1); 
set(gcf,'color','w')



anal_opts.tdc_import.shot_num=find_data_files(anal_opts.tdc_import);

%max_shot_num=max(anal_opts.tdc_import.shot_num);
%anal_opts.tdc_import.shot_num=anal_opts.tdc_import.shot_num(...
%    anal_opts.tdc_import.shot_num>(max_shot_num-anal_opts.history.shots));

data.mcp_tdc=import_mcp_tdc_data(anal_opts.tdc_import);


%% CHECK ATOM NUMBER
sfigure(1)
%create a list of indicies (of the mcp_tdc) that have an ok number of counts
%exclude the very low and then set the thresh based on the sd of the remaining
not_zero_files=data.mcp_tdc.num_counts>1e3; 
num_thresh=mean(data.mcp_tdc.num_counts(not_zero_files))-4*std(data.mcp_tdc.num_counts(not_zero_files));
data.mcp_tdc.num_ok=data.mcp_tdc.num_counts>num_thresh;
 %   (data.mcp_tdc.time_create_write(:,1)'-data.mcp_tdc.time_create_write(1,1))<(anal_opts.max_runtime*60*60);
fprintf('shots number ok %u out of %u \n',sum(data.mcp_tdc.num_ok),numel(data.mcp_tdc.num_ok))

plot((data.mcp_tdc.time_create_write(:,2)-data.mcp_tdc.time_create_write(1,2))/(60*60),data.mcp_tdc.num_counts)
xlabel('time (h)')
ylabel('total counts')
title('num count run trend')


%% Bin pulses

data.mcp_tdc.all_ok=data.mcp_tdc.num_ok;
data.mcp_tdc.al_pulses=bin_al_pulses(anal_opts.atom_laser,data);
%%
anal_opts.osc_fit.adaptive_freq=true; %estimate the starting trap freq 
anal_opts.osc_fit.appr_osc_freq_guess=[52,46.7,40];
anal_opts.osc_fit.plot_fits=true;
anal_opts.osc_fit.plot_err_history=true;
anal_opts.osc_fit.plot_fit_corr=false;

anal_opts.osc_fit.global=anal_opts.global;
data.osc_fit=fit_trap_freq(anal_opts.osc_fit,data);

data.osc_fit.trap_freq_recons=nan*data.osc_fit.ok.did_fits;
mask=data.osc_fit.ok.did_fits;
data.osc_fit.trap_freq_recons(mask)=3*(1/anal_opts.atom_laser.pulsedt)+data.osc_fit.model_coefs(mask,2,1);
data.osc_fit.trap_freq_unc=nan*mask;
data.osc_fit.trap_freq_unc(mask)=data.osc_fit.model_coefs(mask,2,2)';

%%

sfigure(1);
errorbar(data.osc_fit.dld_shot_num,...
    data.osc_fit.trap_freq_recons,data.osc_fit.trap_freq_unc,...
    'kx-','MarkerSize',7,'CapSize',0,'LineWidth',1)
xlabel('Shot Number')
ylabel('Fit Trap Freq')
pause(1e-6)


%% allan


allan_data=[];
allan_data.freq=data.osc_fit.trap_freq_recons;
allan_data.time=data.mcp_tdc.time_create_write(data.osc_fit.dld_shot_idx,2);
mask=~isnan(allan_data.freq) & ~isnan(allan_data.time)';
allan_data.freq=allan_data.freq(mask);
allan_data.time=allan_data.time(mask);
allan_data.time=allan_data.time-min(allan_data.time);
tau_in=2.^(-10:0.05:14)
%%

[retval, s, errorb, tau]=allan(allan_data,tau_in)

%%
windows=[];

fit_out=plot_allan_with_fits(retval, errorb, tau,windows,'Normal Allan Deviation')



%%
[retval, s, errorb, tau]=allan_overlap(allan_data,2.^(-10:0.005:20))
%%
windows=[[57,160];[249,528];[1000,1563];[2091,inf]];
fit_out=plot_allan_with_fits(retval, errorb, tau,windows,'Overlapping Allan Deviation')



%%
[retval, s, errorb, tau]=allan_modified(allan_data,2.^(-10:0.01:20))

%%
windows=[[106,639];[1000,1563];[1885,3092];[4040,inf]];
fit_out=plot_allan_with_fits(retval, errorb, tau,windows,'Modified Allan Deviation')


%%
function fits=plot_allan_with_fits(retval, errorb, tau,windows,dev_label)
font_size = 18;
font_type='cmr10' 
mask_outliers=retval<1;
figure(1)
clf
set(gcf,'color','w')
plot(tau(mask_outliers),retval(mask_outliers),'k-')
hold on
plot(tau(mask_outliers),retval(mask_outliers)+errorb(mask_outliers),'-','Color',[1,1,1]*0.5)
plot(tau(mask_outliers),retval(mask_outliers)-errorb(mask_outliers),'-','Color',[1,1,1]*0.5)
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
pause(1e-6)
xl=xlim;
yl=ylim;
fits=[];
tausamp=10.^linspace(log10(xl(1)),log10(xl(2)),1e3);
labels={dev_label,'+error','-error'};
for ii=1:size(windows,1)
    windows(ii,:)
    mask_fit=tau>windows(ii,1) & tau<windows(ii,2) & mask_outliers;
    fits{ii}=polyfit(log10(tau(mask_fit)),log10(retval(mask_fit)),1);
    plot(tausamp,10.^polyval(fits{ii},log10(tausamp)));
    labels{end+1}=sprintf('Lin. Fit to Section (grad %.2f)',fits{ii}(1))
end
xlim(xl)
ylim(yl)
hold off
title(dev_label,'FontSize',font_size*1.5,'FontName',font_type)
xlabel('\tau (s)','FontSize',font_size,'FontName',font_type)
ylabel('\sigma(\tau) (Hz)','FontSize',font_size,'FontName',font_type)
legend(labels,'FontName',font_type,'FontSize',font_size*0.8,'location','southwest')



end


