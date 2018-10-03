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
%   -osc fit should be weighted properly
%    -freq fits
%    -number fits
%    -number, osc amp corrections
%
% Author: Bryce Henson
% email: Bryce.Henson@live.com
% Last revision:2018-08-19


%close all
%clear all
tic
%%
% BEGIN USER VAR-------------------------------------------------
import_opts.dir='\\amplpc29\Users\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20180814_pulsed_AL_method_tuneout_scan\';
import_opts.file_name='d';
import_opts.force_reimport=false;
import_opts.force_forc=false;
import_opts.dld_xy_rot=0.61;
%Should probably try optimizing these
xlim=[-30e-3, 30e-3];     %tight XY lims to eliminate hot spot from destroying pulse widths
ylim=[-30e-3, 30e-3];
tlim=[.4,1.6];
import_opts.txylim=[tlim;xlim;ylim];

anal_opts.pulsedt=8.000e-3;
anal_opts.t0=0.4178;
anal_opts.start_pulse=1;
anal_opts.pulses=130;
anal_opts.appr_freq_guess=[52,40,40];
anal_opts.pulse_twindow=anal_opts.pulsedt*0.9;
anal_opts.xylim=import_opts.txylim(2:3,:); %set same lims for pulses as import
anal_opts.fall_time=0.417;


histplot.binsx=1000;
histplot.blur=1;
histplot.xlim=[-20,20]*1e-3;
histplot.tlim=[0.86,1.08];
histplot.dimesion=2; %Select coordinate to bin. 1=X, 2=Y.

% END USER VAR-----------------------------------------------------------

%this first section  sets up the struct 'data' which will contain everything you could want incuding the txy data and
%the information from the log
addpath('Colormaps') 
addpath('FileTime_29Jun2011') %used for high precision windows timestamps in import_data
constants
import_opts.shot_num=find_data_files(import_opts);
[data,import_opts]=import_data(import_opts);

%import the log
%adaptively to deal with the 2 different log files that are in the data
clear('log')
log.dir = strcat(import_opts.dir,'log_LabviewMatlab.txt');
fid = fopen(log.dir );
log.cell=textscan(fid,'%s','Delimiter','\n');
log.cell=log.cell{1};
for ii=1:size(log.cell,1)
    if ~isequal(log.cell{ii},'') %catch the empty case
        if contains(log.cell{ii},'measure_probe')
            line_cells=textscan(log.cell{ii},'%f %s %s %s %f %s %u','Delimiter',',');
            log.setpoints(ii)=line_cells{5};
            log.probe_calibration(ii)=false;
            log.iter_nums(ii)=line_cells{7};
        elseif contains(log.cell{ii},'calibrate')
            line_cells=textscan(log.cell{ii},'%f %s %s %s %s %u','Delimiter',',');
            log.setpoints(ii)=NaN;
            log.probe_calibration(ii)=true;
            log.iter_nums(ii)=line_cells{6};
        else %deals with the legacy case (only 20180813_CW_AL_tuneout_scan)
            line_cells=textscan(log.cell{ii},'%f %s %s %s %f %s %u','Delimiter',',');
            log.setpoints(ii)=line_cells{5};
            log.probe_calibration(ii)=false;
            log.iter_nums(ii)=line_cells{7};
        end
        log.posix_times(ii)=line_cells{1};
        log.iso_times{ii}=line_cells{2};
    end
end
data.probe.setpoint=log.setpoints*1e6; %convert to hz
data.probe.time=log.posix_times;
data.probe.shot_num=log.iter_nums;
data.probe.calibration=log.probe_calibration;
%can check that the times look ok
% plot(data.mcp_tdc.write_time-data.probe.time)

%create a list of indicies that have an ok number of counts
%exclude the very low and then set the thresh based on the sd of the remaining
not_zero_files=data.mcp_tdc.num_counts>1e2; 
num_thresh=mean(data.mcp_tdc.num_counts(not_zero_files))-3*std(data.mcp_tdc.num_counts(not_zero_files));
data.mcp_tdc.num_ok=data.mcp_tdc.num_counts>num_thresh;


%check that the shot numbers align
if ~isequal(data.probe.shot_num,data.mcp_tdc.shot_num),error('shot numbers do not align!'),end

%create a list of shots that belong to a particular probe beam set pt
unique_freqs=unique(data.probe.setpoint);
unique_freqs=unique_freqs(~isnan(unique_freqs));
for ii=1:numel(unique_freqs)
    freq=unique_freqs(ii);
    %get the itteration numbers that had wavelength equaling freq
    matching_itt_num=data.probe.shot_num(freq==data.probe.setpoint);
    data.freq_sorted.freq(ii)=freq;
    data.freq_sorted.indexes{ii}=matching_itt_num;
    %I think this is cleaner than copying data arround a lot, however if there is a lot of reporcessing going on the
    %line below can just write it out
    %data.freq_sorted.txy={data.mcp_tdc.counts_txy{matching_itt_num}};
end
%usage example
%to get all the files that are at freq =  data.freq_sorted.freq(5)
%{data.mcp_tdc.counts_txy{data.freq_sorted.indexes{5}}}
%gives a cell array each of which contains matricies of a shots txy data
%vertcat(data.mcp_tdc.counts_txy{data.freq_sorted.indexes{5}})
%gives all the data in a given wavelength


%% binning the pulses
%first step is to bin up each pulse of the atom laser

tic
data.mcp_tdc.al_pulses=[];
data.mcp_tdc.al_pulses.pulsedt=anal_opts.pulsedt;
data.mcp_tdc.al_pulses.window=nan(anal_opts.pulses,3,2); %initalize
data.mcp_tdc.al_pulses.num_counts=nan(size(data.mcp_tdc.counts_txy,2),anal_opts.pulses);
  
plots=false;
pulse_win_txy_arr=[];
fprintf('Binning pulses in %04i shots\n%04i\n',size(data.mcp_tdc.counts_txy,2),0)
for shot=1:size(data.mcp_tdc.counts_txy,2)
        for pulse=1:anal_opts.pulses
            %set up time window centered arround t0
            trange=anal_opts.t0+anal_opts.pulsedt*(anal_opts.start_pulse+pulse-2)+anal_opts.pulse_twindow*[-0.5,0.5];
            pulse_win_txy=[trange;anal_opts.xylim]; 
            counts_pulse=masktxy(data.mcp_tdc.counts_txy{shot},pulse_win_txy);
            if plots
                figure(1)
                set(gcf,'Color',[1 1 1]);
                subplot(3,1,1)
                hist(counts_pulse(:,1),100)
                xlabel('t')
                title('full')
                subplot(3,1,2)
                hist(counts_pulse(:,2),100)
                xlabel('x')
                title('sub')
                subplot(3,1,3)
                hist(counts_pulse(:,3),100)
                xlabel('y')
                title('sub')
                pause(0.01)
            end
            %this isnt really needed most of the time
            data.mcp_tdc.al_pulses.counts_txy{shot,pulse}=counts_pulse;
            if shot==1 
                data.mcp_tdc.al_pulses.window(pulse,:,:)=pulse_win_txy; %only need to store this on first shot
                data.mcp_tdc.al_pulses.time(pulse,:)=(pulse+anal_opts.start_pulse-2)*anal_opts.pulsedt;
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
  fprintf('\b\b\b\b%04u',shot)     
%to set the pulse t0 right it can be handy to uncomment the next line
%fprintf('\nmean time %3.5f            \n ',mean(anal_out.mean_pos_window(file,:,3)-anal_out.mean_pos_window(file,:,1)))
%anal_out=single_seg_fig(data,import_opts,anal_out,anal_opts,histplot,4,file);
end%shots
fprintf('\ndone Binning\n')

toc


%% Fitting the Trap freq
%using the binned data we fit the trap freq to each shot

% deltphase=0.05*2*pi;
% fprintf('xy angle %2.5f° \n',asin(fitparam{5,1})*(180/pi))
% fprintf('xz angle %2.5f° \n',asin(fitparam{6,1})*(180/pi))
% fprintf('time untill phase unc reaches %2.3f°  %2.3fcyc= %2.3f \n',...
%     deltphase*(180/pi),deltphase/(2*pi),(deltphase-fitparam{3,2})/fitparam{2,2})




plots=false;
anal_opts.appr_freq_guess=[52,47.6,40];

velocity=const.g0*anal_opts.fall_time;
warning('off','stats:nlinfit:ModelConstantWRTParam');
warning('off','MATLAB:rankDeficientMatrix');
data.osc_fit=[]; %clear the output struct

good_num_idx=1:numel(data.mcp_tdc.num_ok);
good_num_idx=good_num_idx(data.mcp_tdc.num_ok);
fprintf('Fitting oscillations in %04i shots\n%04i\n',size(good_num_idx,2),0)
for ii=1:size(good_num_idx,2)
shot=good_num_idx(ii);
data.osc_fit.shot(ii)=shot;
%construct a more convinent temp variable txyz_tmp wich is the position in mm for use in the fit
x_tmp=1e3*squeeze(data.mcp_tdc.al_pulses.pos_stat(shot,:,2));
x_tmp=x_tmp-mean(x_tmp);
y_tmp=1e3*squeeze(data.mcp_tdc.al_pulses.pos_stat(shot,:,3));
y_tmp=y_tmp-mean(y_tmp);
z_tmp=data.mcp_tdc.al_pulses.time-squeeze(data.mcp_tdc.al_pulses.pos_stat(shot,:,1))';
z_tmp=z_tmp'*velocity*1e3;
z_tmp=z_tmp-mean(z_tmp);
txyz_tmp=[data.mcp_tdc.al_pulses.time';x_tmp;y_tmp;z_tmp];
sqrtn=sqrt(data.mcp_tdc.al_pulses.num_counts(shot,:)); %find the statistical uncert in a single shot
sqrtn=sqrtn;
xerr=1e3*squeeze(data.mcp_tdc.al_pulses.pos_stat(shot,:,5))./sqrtn;
yerr=1e3*squeeze(data.mcp_tdc.al_pulses.pos_stat(shot,:,6))./sqrtn;
zerr=1e3*squeeze(data.mcp_tdc.al_pulses.pos_stat(shot,:,4))*velocity./sqrtn;
xyzerr_tmp=[xerr;yerr;zerr];

%modelfun = @(b,x) b(1)*sin(b(2)*x(:,1)*pi+b(3)*pi)+b(4)+b(5)*x(:,1)+x(:,2)*b(6)+x(:,3)*b(7);
%beta0=[15, 30, -0.5, -2, 0, 0, 0];
%cof_names={'amp','freq','phase','offset','grad','ycpl','zcpl'}

modelfun = @(b,x) exp(-x(:,1).*max(0,b(7))).*b(1).*sin(b(2)*x(:,1)*pi*2+b(3)*pi*2)+b(4)+b(8)*x(:,1)+b(5)*x(:,2)*0+b(6)*x(:,3)*0;
beta0=[5, anal_opts.appr_freq_guess(histplot.dimesion), 1, 1,0,0,0.5,0.1];
cof_names={'amp','freq','phase','offset','ycpl','zcpl','damp','grad'};

opt = statset('TolFun',1e-10,...
    'TolX',1e-10,...
    'MaxIter',1e4,...
    'UseParallel',1);

%select the aproapriate values to go in the response variable
idx=1:4;
idx(histplot.dimesion+1)=[];
predictor=txyz_tmp(idx,:)';
%predictor=[tvalues,xvalues,zvalues];
fitobject=fitnlm(predictor,txyz_tmp(histplot.dimesion+1,:)',modelfun,beta0,'options',opt,...
    'CoefficientNames',cof_names);
data.osc_fit.model{shot}=fitobject;
fitparam=fitobject.Coefficients;
data.osc_fit.model_coefs(ii,:,:)=[a.Coefficients.Estimate,a.Coefficients.SE];
data.osc_fit.fit_rmse(ii)=a.RMSE;
%limiting frequnecy prediction from http://adsabs.harvard.edu/full/1999DSSN...13...28M
meanwidth=sqrt(mean(squeeze(data.mcp_tdc.al_pulses.pos_stat(shot,:,5)).^2))*1e3;
frequnclim=sqrt(6/sum(data.mcp_tdc.al_pulses.num_counts(shot,:)))*...
    (1/(pi*range(data.mcp_tdc.al_pulses.time)))*...
    (meanwidth/fitparam{2,1});
%fprintf('sampling limit %2.3g Hz, fint unc %2.3g Hz, ratio %2.3g \n',[frequnclim,fitparam{2,2},fitparam{2,2}/frequnclim])
data.osc_fit.fit_sample_limit{shot}=[frequnclim,fitparam{2,2},fitparam{2,2}/frequnclim];

tplotvalues=linspace(min(data.mcp_tdc.al_pulses.time),max(data.mcp_tdc.al_pulses.time),1e5)';
predictorplot=[tplotvalues,...
               interp1(predictor(:,1),predictor(:,2),tplotvalues),...
               interp1(predictor(:,1),predictor(:,3),tplotvalues)];
[prediction,ci]=predict(fitobject,predictorplot);
%set(gcf, 'Units', 'Pixels', 'OuterPosition', [100, 100, 500, 300])
%set(gca,'ylim',[-20,20])
%yl=ylim;
%yticks(linspace(yl(1),yl(2),6))
%yticks([-8 -4 0 4 8])
pause(0.01)


if plots
    sfigure(2)
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
    pause(0.1)

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
    saveas(gca,sprintf('.\\out\\fit_shot%04u.png',shot))
end
fprintf('\b\b\b\b%04u',ii)
end
fprintf('\ndone Fitting\n')



%%
%now extract what we care about from the fit objects
for ii=1:size(good_num_idx,2)
shot=good_num_idx(ii);
data.osc_fit.shot(ii)=shot;
a=data.osc_fit.model{shot};
data.osc_fit.model_coefs(ii,:,:)=[a.Coefficients.Estimate,a.Coefficients.SE];
data.osc_fit.fit_rmse(ii)=a.RMSE;
end
%%
%look for failed fits
subplot(2,1,1)
hist(data.osc_fit.fit_rmse)
subplot(2,1,2)
plot(data.osc_fit.fit_rmse)

data.osc_fit.good_fit=data.osc_fit.fit_rmse< mean(data.osc_fit.fit_rmse)+2*std(data.osc_fit.fit_rmse);
data.osc_fit.good_fit=data.osc_fit.good_fit & ...
    abs(data.osc_fit.model_coefs(:,2,1)'-mean(data.osc_fit.model_coefs(:,2,1)))<...
    5;
%%
plot(data.probe.setpoint(data.osc_fit.shot(data.osc_fit.good_fit)),data.osc_fit.model_coefs(data.osc_fit.good_fit,2,1),'x')


%%
%create a mode of the underlying frequency
figure(3)
clf
set(gcf,'color','w')
idx=data.probe.calibration(data.osc_fit.shot);
x_tmp=data.mcp_tdc.write_time(idx);
time_start=x_tmp(1);
x_tmp=x_tmp-time_start;
y_tmp=3*(1/anal_opts.pulsedt)+data.osc_fit.model_coefs(idx,2,1);
freq_drift_model=fit(x_tmp', y_tmp,'poly9');
x_samp=linspace(min(x_tmp),max(x_tmp),1e3);
hour_in_s=60*60;
plot(x_tmp/hour_in_s,y_tmp)
hold on
plot(x_samp/hour_in_s,freq_drift_model(x_samp))
xlabel('experiment time (h)')
ylabel('no probe trap freq')
legend('data','polynomial model')


%%
figure(4)
clf
subplot(2,1,1)
probe_freq=data.probe.setpoint(data.osc_fit.shot(data.osc_fit.good_fit));
corrected_delta=data.osc_fit.model_coefs(data.osc_fit.good_fit,2,1)-...
    freq_drift_model(data.mcp_tdc.write_time(data.osc_fit.shot(data.osc_fit.good_fit))-time_start);
plot(probe_freq,corrected_delta,'x')
subplot(2,1,2)
square_diff= (3*(1/anal_opts.pulsedt)+data.osc_fit.model_coefs(data.osc_fit.good_fit,2,1)).^2-...
    (freq_drift_model(data.mcp_tdc.write_time(data.osc_fit.shot(data.osc_fit.good_fit))-time_start)).^2;
plot(probe_freq,square_diff,'x')


%%

%%
%now we do some fitting
%thresholds for CI
%sd         CI
%1          0.3174
%2          0.05
%3          2.699e-03
ci_size_disp=0.3174;%one sd %confidence interval to display
ci_size_cut_outliers=0.05; %2.699e-03= 3sigma %confidence interval to kill outliers

fprintf('Calculating Fits')
%select the data in some freq range and that has an ok number

%set up the data input for the fit
xdat=probe_freq(~isnan(probe_freq))';
ydat=square_diff(~isnan(probe_freq))';
if exist('new_to_freq_val','var') %if the TO has been calculated before use that as the center
    freq_offset=new_to_freq_val;
else
    freq_offset=nanmean(xdat); %otherwise use the center of the range
end
xdat=xdat-freq_offset; %fits work better when thery are scaled reasonably
xdat=xdat*1e-9;

modelfun = @(b,x) b(1)+ b(2).*x; %simple linear model
opts = statset('nlinfit');
opts.RobustWgtFun = 'welsch' ; %a bit of robust fitting
opts.Tune = 1;
beta0 = [1e-5,1e-2]; %intial guesses
mdl = fitnlm(xdat,ydat,modelfun,beta0,'Options',opts)
xsamp=linspace(min(xdat),max(xdat),1e3)'; %sample for the model curve
[ysamp,yci]=predict(mdl,xsamp,'Prediction','observation','Alpha',ci_size_disp); %note the observation CI
%now plot the data and the model together
figure(6);
set(gcf,'color','w')
subplot(1,2,1)
plot(xsamp,ysamp,'k-')
hold on
plot(xsamp,yci,'r-')
plot(xdat,ydat,'bx')
xlabel(sprintf('probe beam set freq - %.3f(GHz)',freq_offset*1e-9))
ylabel('Response (Hz^2)')
title('Good Data')
first_plot_lims=[get(gca,'xlim');get(gca,'ylim')];
%now make a new prediction with the model but with the CI to cut out outliers
[~,yci_cull_lim]=predict(mdl,xdat,'Prediction','observation','Alpha',ci_size_cut_outliers);
color_idx=ydat>yci_cull_lim(:,1)' & ydat<yci_cull_lim(:,2)';
%color the ones that will be removed
plot(xdat(~color_idx),ydat(~color_idx),'r.','markersize',15)
hold off
xdat_culled=xdat(color_idx);
ydat_culled=ydat(color_idx);


opts = statset('nlinfit');
%opts.RobustWgtFun = 'welsch' ;
opts.Tune = 1;
beta0 = [1e-5,1e-2];
mdl_culled = fitnlm(xdat_culled,ydat_culled,modelfun,beta0,'Options',opts)
xsamp_culled=linspace(min(xdat_culled),max(xdat_culled),1e3)';
[ysamp_culled,yci_culled]=predict(mdl_culled,xsamp_culled,'Alpha',0.2); %'Prediction','observation'
%now plot the remaining data along with the fit model and the model CI
subplot(1,2,2)
plot(xsamp_culled,ysamp_culled,'k-')
hold on
plot(xsamp_culled,yci_culled,'r-')
plot(xdat_culled,ydat_culled,'bx')
hold off
xlabel(sprintf('probe beam set freq - %.3f (GHz)',freq_offset*1e-9))
ylabel('Response (Hz^2)')
title('Fit Outliers Removed')
set(gca,'xlim',first_plot_lims(1,:))
set(gca,'ylim',first_plot_lims(2,:))

%normalize by the CI at the TO
figure(7);
clf
set(gcf,'color','w')
cross_xval=-mdl.Coefficients.Estimate(1)/mdl.Coefficients.Estimate(2);
[cross_yval,cross_yci]=predict(mdl,cross_xval,'Prediction','observation','Alpha',ci_size_disp);
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
%inverse scaled gradient to give the single shot uncert
single_shot_uncert=abs(1/(mdl.Coefficients.Estimate(2)/cross_yci));
fprintf('single shot uncert detuning @1SD %.1f MHz, %.2f fm\n',single_shot_uncert*1e3,...
    1e9*single_shot_uncert*const.c/((new_to_freq_val*2)^2)*10^15)
%predicted uncert using this /sqrt(n)
fprintf('predicted stat. uncert %.1f MHz, %.2f fm\n',single_shot_uncert/sqrt(sum(data.mcp_tdc.num_ok))*1e3,...
    1e9*single_shot_uncert/sqrt(sum(data.mcp_tdc.num_ok))*const.c/((new_to_freq_val*2)^2)*10^15)

fprintf('\nfit results\n')
fprintf('median damping time %.2f\n',median(1./data.osc_fit.model_coefs(data.osc_fit.good_fit,7,1)))
%calculate some statistics and convert the model parameter into zero crossing and error therin
old_to_wav=413.0938e-9;
new_to_freq_val=-1e9*mdl_culled.Coefficients.Estimate(1)/mdl_culled.Coefficients.Estimate(2)+freq_offset;
new_to_freq_unc=1e9*abs((mdl_culled.Coefficients.Estimate(1)/mdl_culled.Coefficients.Estimate(2)))*...
    sqrt((mdl_culled.Coefficients.SE(1)/mdl_culled.Coefficients.Estimate(1))^2+...
    (mdl_culled.Coefficients.SE(2)/mdl_culled.Coefficients.Estimate(2))^2);
new_to_wav_val=const.c/(new_to_freq_val*2);
new_to_wav_unc=new_to_freq_unc*const.c/((new_to_freq_val*2)^2);
fprintf('TO freq        %.2f±%.2f MHz\n',new_to_freq_val*1e-6,new_to_freq_unc*1e-6)
fprintf('TO wavelength  %.6f±%f nm \n',new_to_wav_val*1e9,new_to_wav_unc*1e9)
fprintf('diff from TOV1 %e±%e nm \n',(new_to_wav_val-old_to_wav)*1e9,new_to_wav_unc*1e9)
fprintf('number of probe files        %u \n',sum(~isnan(probe_freq)))
fprintf('number of calibration files  %u \n',sum(isnan(probe_freq)))
fprintf('total used                   %u \n',size(probe_freq,2))
fprintf('files with enough number     %u\n',sum(data.mcp_tdc.num_ok))

%% damping results
%plot out what the distibution over damping times is
figure(7);
set(gcf,'color','w')
histogram(1./data.osc_fit.model_coefs(data.osc_fit.good_fit,7,1),linspace(0,10,1e2))

toc