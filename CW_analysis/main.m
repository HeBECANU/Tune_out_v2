%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2018-08-14
% Author: J Ross, Bryce M Henson (bryce.henson@live.com)
% A script to analyse the data we collected from a CW atom laser whose output flux is modulated by application of a probe beam near the tuneout wavelength.
% The script:
%   * Imports all the data files into a data structure
%   * Imports the log file, which indexes all the data as a) data run,
%       including wavelength, or b) calibration run to check the control
%       signal
%   * extracts the response at each set point using a lock-in method
%   * Produces a plot of the measured response vs probe
%       wavelength, identifying the tuneout wavelength and giving a
%       stastistical uncertainty.
%the data structure
%   first level is instrument or anal method
%   try and pass things between the modules of code only using the 'data' struct

%Known BUGS/ possible Improvements--------------------------------------------------------
%2018-08-16
%could quite happily remove the FFT stuff
%look at how the CW al TOF profile changes throughut the run
%   *can dramaticaly change the sensitivity which will make the spread in the TO plot much hihger

%Fixed !
%2018-08-16
%putting all the files in the apropriate main directory
%user var section seperate from code
%new import code used (faster)
%fixed histgraming error that was off by 1 bin
%fixed steady progression of wavelength assumption
%fixed fft_tx
%made data structure for freq indexes
%made the data struct more hierarchical so it is better self documenting eg
%   data.mcp_tdc
%   data.probe


%close all
%clear alls
tic
%%
% BEGIN USER VAR-------------------------------------------------

%import_opts.dir='\\amplpc29\Users\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20180813_CW_AL_tuneout_scan\';
%import_opts.dir='\\amplpc29\Users\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20180814_CW_AL_method_polz_dependence\no_wp\';
%import_opts.dir='\\amplpc29\Users\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20180814_CW_AL_method_polz_dependence\wp_289deg_good_align';
%import_opts.dir='\\amplpc29\Users\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20180814_CW_AL_method_polz_dependence\wp_289deg_no_align';
%import_opts.dir='\\amplpc29\Users\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20180814_CW_AL_method_polz_dependence\wp_312deg_bad_align\';
%import_opts.dir='\\amplpc29\Users\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20180814_CW_AL_method_polz_dependence\wp_312deg_good_align\';
%import_opts.dir='\\amplpc29\Users\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20180814_CW_AL_method_polz_dependence\wp_223deg_good_align';
%import_opts.dir='\\amplpc29\Users\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20180814_CW_AL_method_polz_dependence\wp_267deg_good_align';
import_opts.dir='\\amplpc29\Users\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\';
import_opts.file_name='d';
import_opts.force_reimport=false;
import_opts.force_forc=false;
import_opts.dld_xy_rot=0.61;
%Should probably try optimizing these
xlim=[-20e-3, 20e-3];     %tight XY lims to eliminate hot spot from destroying pulse widths
ylim=[-20e-3, 20e-3];
tlim=[.4,1.4];
import_opts.txylim=[tlim;xlim;ylim];

anal_opts.hist.tbin = 5e-5;
anal_opts.hist.tlim = import_opts.txylim(1,:);
anal_opts.hist.tsmooth=0e-6;%2e-4;
anal_opts.mod_freq=427;
anal_opts.do_fft=false;
anal_opts.plot_fft=false;
% END USER VAR-----------------------------------------------------------
%this first section  sets up the struct 'data' which will contain everything you could want incuding the txy data and
%the information from the log
addpath('Colormaps') 
addpath('FileTime_29Jun2011') %used for high precision windows timestamps in import_data
constants
import_opts.shot_num=find_data_files(import_opts);
import_opts.shot_num=import_opts.shot_num(import_opts.shot_num<400) %HACK
[data,import_opts]=import_data(import_opts);
if ~exist([import_opts.dir,'out\'],'dir')
    mkdir([import_opts.dir,'out\'])
end

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

%create a list of indicies that have an ok number of counts
%exclude the very low and then set the thresh based on the sd of the remaining
not_zero_files=data.mcp_tdc.num_counts>1e2; 
num_thresh=mean(data.mcp_tdc.num_counts(not_zero_files))-3*std(data.mcp_tdc.num_counts(not_zero_files));
data.mcp_tdc.num_ok=data.mcp_tdc.num_counts>num_thresh;


aligned_idx=min([size(data.probe.shot_num,2),size(data.mcp_tdc.shot_num,2),size(data.mcp_tdc.shot_num,2)]);
data.probe.shot_num=data.probe.shot_num(1:aligned_idx);
data.mcp_tdc.shot_num=data.mcp_tdc.shot_num(1:aligned_idx);
data.probe.setpoint=data.probe.setpoint(1:aligned_idx);
data.probe.calibration=data.probe.calibration(1:aligned_idx);
%check that the shot numbers align
if ~isequal(data.probe.shot_num,data.mcp_tdc.shot_num),error('shot numbers do not align'),end

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


%%

%calculate the histograms and save the counts in the  data.mcp_tdc.hist struct
 data.mcp_tdc.hist=[];%clear this part of the stuct
 data.mcp_tdc.hist.opts=anal_opts.hist;
%now set up and calculate the histogram
%general problem with histograming is when the time window is not an exact integer number of desired bin size wide
%the least bad solution is to round up the number of bins and have a half full bin at the end
num_shots = length(import_opts.shot_num);
num_bins = ceil(diff(anal_opts.hist.tlim)/anal_opts.hist.tbin); %round up 
hist_bin_edges = linspace(import_opts.txylim(1,1),import_opts.txylim(1,1)+num_bins*anal_opts.hist.tbin,num_bins);
%previously this assumed that the data was constant 
%hist_bin_centers = hist_bin_edges(1:end-1) + 0.5*t_stepsize;
 data.mcp_tdc.hist.bin_centers = 0.5*(hist_bin_edges(1:end-1) + hist_bin_edges(2:end));

 data.mcp_tdc.hist.counts_raw=nan(size(data.mcp_tdc.counts_txy,2),size( data.mcp_tdc.hist.bin_centers,2)); %prealocate
if anal_opts.hist.tsmooth~=0
     data.mcp_tdc.hist.counts_smooth= data.mcp_tdc.hist.counts_raw;
end
fprintf('Calculating %04i histograms \n%04i\n',size(data.mcp_tdc.counts_txy,2),0)
data.lock_in.quad_raw=nan(size(data.mcp_tdc.counts_txy,2),2);%initalize with nan so failure is obvious
for ii=1:size(data.mcp_tdc.counts_txy,2)
    fprintf('\b\b\b\b%04u',ii)
     data.mcp_tdc.hist.counts_raw(ii,:)=histcounts(data.mcp_tdc.counts_txy{ii}(:,1),hist_bin_edges)/anal_opts.hist.tbin;
    if anal_opts.hist.tsmooth~=0
        data.mcp_tdc.hist.counts_smooth(ii,:)=gaussfilt(hist_bin_centers, data.mcp_tdc.hist.counts_raw(ii,:),anal_opts.hist.tsmooth);
    end
end
fprintf('\nHistograms Done\n')



%ok now it looks like we are ready to try and see some results...
%%
%the first way i tried doing the detection of the signal was to use the FFT
%this is nice in that it can show you interesting features in the spectrum
%I then took the phasor from the fft and then rotated it so the signal was entirely in the real part
%this real part is then the signal
%the advantage of this over previous method that took the magnitude of the pasor 
%   *is that the average over noise is ZERO
%   *there is none of this fuckery where you multiply by -ve one
%   *there is (at best) a sqrt(2) improvement in noise because now the out of phase component does not add
%note because i moved on to the lock in method i have not implemented the auto phase alignment that the lock in method
%has
if anal_opts.do_fft
    fprintf('taking fft\n%04i\n',0)
    for ii=1:size(data.mcp_tdc.counts_txy,2)
        if anal_opts.hist.tsmooth==0
           hist_counts= data.mcp_tdc.hist.counts_raw(ii,:);
        else
           hist_counts= data.mcp_tdc.hist.counts_smooth(ii,:);
        end
        fftout=fft_tx( data.mcp_tdc.hist.bin_centers,hist_counts,10);%use a padded fft for smaller freq bins
        [~,index]=min(abs(fftout(1,:)-anal_opts.mod_freq));
        %sample_freq=fftout(1,index)
        data.mod_amp_fft(ii)=fftout(2,index);
        if anal_opts.plot_fft
            sfigure(1);
            clf
            subplot(1,2,1)
            plot( data.mcp_tdc.hist.bin_centers,hist_counts) 
            subplot(1,2,2)
            plot(fftout(1,:),abs(fftout(2,:)))
            set(gca,'xlim',[400,450]) 
            pause(0.01)
        end
        fprintf('\b\b\b\b%04u',ii)
    end
    fprintf('\ndone fft\n')

    %%
    %plot the magnitude and phase of the phasor
    %you should see a v shape in the magnitude and a step in the phase
    sfigure(2);
    set(gcf,'color','w')
    subplot(1,2,1)
    plot(data.probe.setpoint(data.mcp_tdc.num_ok),abs(data.mod_amp_fft(data.mcp_tdc.num_ok)),'xk')
    xlabel('pobe beam set freq')
    ylabel('response magnitude')
    subplot(1,2,2)
    plot(data.probe.setpoint(data.mcp_tdc.num_ok),angle(data.mod_amp_fft(data.mcp_tdc.num_ok)),'xk')
    xlabel('pobe beam set freq')
    ylabel('response angle') 

    %%
    %lets try rotating it so the signal is in the real axis
    %find the mean angle at the start
    %just need some data that is one side of the TO to estimate the roation angle
    indicies=data.probe.setpoint<3.62855e14 & data.mcp_tdc.num_ok;
    mean_angle=mean(angle(data.mod_amp_fft(indicies)));
    theta = mean_angle-0.45*pi+pi/2; % to rotate 90 counterclockwise
    %i do this rotation by conversion into vectors wich is equivelent and slower to the method below
    vector_rep=[real(data.mod_amp_fft);imag(data.mod_amp_fft)];
    %dot(vector_rep(:,1),in_phase_unit)
    R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    %create unit vectors
    in_phase_unit = R*([1 0]');
    out_phase_unit = R*([0 1]');

    in_phase_comp=dot(repmat(in_phase_unit,[1,size(vector_rep,2)]),vector_rep,1);
    out_phase_comp=dot(repmat(out_phase_unit,[1,size(vector_rep,2)]),vector_rep,1);

    figure(3);
    clf
    set(gcf,'color','w')
    subplot(1,2,1)
    plot(data.probe.setpoint(data.mcp_tdc.num_ok),in_phase_comp(data.mcp_tdc.num_ok),'x')
    xlabel('pobe beam set freq')
    ylabel('In phase response')
    subplot(1,2,2)
    plot(data.probe.setpoint(data.mcp_tdc.num_ok),out_phase_comp(data.mcp_tdc.num_ok),'x')
    xlabel('pobe beam set freq')
    ylabel('Out phase response')



    %%
    %you can just rotate about the complex plane by simple multipication wich is heaps easier
    cpx_rep=exp(-theta*1j).*data.mod_amp_fft;
    figure(3);
    clf;
    set(gcf,'color','w')
    subplot(1,2,1)
    plot(data.probe.setpoint(data.mcp_tdc.num_ok),real(cpx_rep(data.mcp_tdc.num_ok)),'x')
    xlabel('pobe beam set freq')
    ylabel('In phase response')
    subplot(1,2,2)
    plot(data.probe.setpoint(data.mcp_tdc.num_ok),imag(cpx_rep(data.mcp_tdc.num_ok)),'x')
    xlabel('pobe beam set freq')
    ylabel('Out phase response')
    mean(angle(cpx_rep(indicies)))

end
%%
%compute using lock in
%the lock in method multiplies by sin and cos at the freq of interest
%it is a LOT faster than the FFT in principle(O(n) instead of O(n log(n))
%by precomputing the histograms and sin/cos arrays and just using element wise multiplication everything is insanely
%quick
fprintf('calculating Lock In...')
precompute_trig=[sin(2*pi*anal_opts.mod_freq* data.mcp_tdc.hist.bin_centers);...
                cos(2*pi*anal_opts.mod_freq* data.mcp_tdc.hist.bin_centers)];
data.lock_in.freq=anal_opts.mod_freq;  
%make this have dimension [shots,bins,2]
precompute_trig=repmat(permute(precompute_trig, [3 2 1]),[size( data.mcp_tdc.hist.counts_raw,1),1,1]);
%then repeat the data twice in the 3rd dim

if anal_opts.hist.tsmooth~=0
    lock_in_data= data.mcp_tdc.hist.counts_smooth;
else
     lock_in_data= data.mcp_tdc.hist.counts_raw;
end
lock_in_data=repmat(lock_in_data,[1,1,2]);

lock_in_data=lock_in_data.*precompute_trig;
lock_in_data=permute(sum(lock_in_data,2),[1,3,2]);
lock_in_data=lock_in_data/(range( data.mcp_tdc.hist.bin_centers)*numel( data.mcp_tdc.hist.bin_centers));
data.lock_in.quad_raw=lock_in_data;
fprintf('Done\n')

%%
%agian we caluclate the average angle on one side of the tune out where the signal is small
indicies=data.probe.setpoint<3.62855e14;
mean_li_angle=mean(angle(data.lock_in.quad_raw (indicies,1)+data.lock_in.quad_raw (indicies,2)*1j));

%this is only an ok guess, to make it perfect we change the rotation angle while calulating the
%out of phase signal standard deviation
%this gives a very clear minima
dtheta=linspace(-0.5,0.5,1e4);
out_phase_std=nan*dtheta; %initalize output array with nans 
plots=false; %you can watch if you like
fprintf('Calculating Phase Alingment...')
for ii=1:numel(dtheta)
theta =mean_li_angle+dtheta(ii)*pi;%0.045-pi/2; % to rotate 90 counterclockwise
li_data=data.lock_in.quad_raw ;
M=exp(-theta*1j).*(li_data(:,1)+li_data(:,2)*1j);
li_data = [real(M), imag(M)];
if plots
    sfigure(4);
    clf
    set(gcf,'color','w')
    subplot(1,2,1)
    plot(data.probe.setpoint(data.mcp_tdc.num_ok),li_data(data.mcp_tdc.num_ok,1),'x')
    xlabel('pobe beam set freq')
    ylabel('In phase response')
    title(sprintf('%.2f pi',dtheta(ii)))
    subplot(1,2,2)
    plot(data.probe.setpoint(data.mcp_tdc.num_ok),li_data(data.mcp_tdc.num_ok,2),'x')
    xlabel('pobe beam set freq')
    ylabel('Out phase response')
    pause(0.001)
end
out_phase_std(ii)=std(li_data(data.mcp_tdc.num_ok,2));
end

figure(44)
set(gcf,'color','w')
plot(mean_li_angle+dtheta*pi,out_phase_std)
xlabel('rotation')
ylabel('Out phase signal std')
saveas(gcf,[import_opts.dir,'out\phase_align.png']) 
[~,min_idx]=min(out_phase_std);
optimized_theta=mean_li_angle+dtheta(min_idx)*pi;
fprintf('Done\n')

theta =optimized_theta;%0.045-pi/2; % to rotate 90 counterclockwise
M=exp(-theta*1j).*(data.lock_in.quad_raw(:,1)+data.lock_in.quad_raw(:,2)*1j);
data.lock_in.quad_rot=[real(M), imag(M)];
data.lock_in.rot_angle=optimized_theta;

%%
figure(4);
clf
set(gcf,'color','w')
subplot(1,2,1)
plot(data.probe.setpoint(data.mcp_tdc.num_ok),data.lock_in.quad_rot(data.mcp_tdc.num_ok,1),'x')
xlabel('pobe beam set freq')
ylabel('In phase response')
title(sprintf('%.2f pi',dtheta))
subplot(1,2,2)
plot(data.probe.setpoint(data.mcp_tdc.num_ok),data.lock_in.quad_rot(data.mcp_tdc.num_ok,2),'x')
xlabel('pobe beam set freq')
ylabel('Out phase response')
pause(0.1)
saveas(gcf,[import_opts.dir,'out\aligned_response.png']) 

%%
%calibration data
%check that the mean is zero and the std is small
fprintf('caibration data signal mean %.2e (in) %.2e(out)\n',...
    mean(data.lock_in.quad_rot(data.mcp_tdc.num_ok&data.probe.calibration,:)))
%and then compare the std of the calibration signal to the probe data
fprintf('                       std  %.2e (in) %.2e(out)\n',...
    std(data.lock_in.quad_rot(data.mcp_tdc.num_ok&data.probe.calibration,:)))
fprintf('probe data signal      std  %.2e (in) %.2e(out)\n',...
    mean(cellfun(@(idx) std(data.lock_in.quad_rot(idx,1)),data.freq_sorted.indexes)),...
    mean(cellfun(@(idx) std(data.lock_in.quad_rot(idx,2)),data.freq_sorted.indexes)))
%ok so it looks like the modulation increases the variability above the noise floor for the in phase signal

%%
%now that i have a good method of finding the signal
%try and cut out those nasty outliers
figure(5);
plot(data.lock_in.quad_rot(data.mcp_tdc.num_ok,1))
%ok so there is some kind of structure in the ime data but not clear enough to get rid of
%maybe just use a robust fit to the data

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
indicies=data.probe.setpoint>0 & data.mcp_tdc.num_ok;%3.6286e14;
%set up the data input for the fit
xdat=data.probe.setpoint(indicies)';
if exist('new_to_freq_val','var') %if the TO has been calculated before use that as the center
    freq_offset=new_to_freq_val;
else
    freq_offset=mean(xdat); %otherwise use the center of the range
end
xdat=xdat-freq_offset; %fits work better when thery are scaled reasonably
xdat=xdat*1e-9;
ydat=data.lock_in.quad_rot(indicies,1)';
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
ylabel('Response (Hz)')
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
ylabel('Response (Hz)')
title('Fit Outliers Removed')
set(gca,'xlim',first_plot_lims(1,:))
set(gca,'ylim',first_plot_lims(2,:))
saveas(gcf,[import_opts.dir,'out\TO_fit.png']) 

%normalize by the CI at the TO
figure(7);
set(gcf,'color','w')
cross_xval=-mdl.Coefficients.Estimate(1)/mdl.Coefficients.Estimate(2);
[cross_yval,cross_yci]=predict(mdl,cross_xval,'Prediction','observation','Alpha',ci_size_disp);
if abs(cross_yval)>1e-5, error('not crossing zero here') ,end
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
fprintf('single shot uncert @1SD %.1e GHz\n',single_shot_uncert)
%predicted uncert using this /sqrt(n)
%single_shot_uncert/sqrt(sum(data.mcp_tdc.num_ok))
saveas(gcf,[import_opts.dir,'out\TO_sens.png']) 

fprintf('\nfit results\n')
%calculate some statistics and convert the model parameter into zero crossing and error therin
old_to_wav=413.0938e-9;
new_to_freq_val=-1e9*mdl.Coefficients.Estimate(1)/mdl.Coefficients.Estimate(2)+freq_offset;
new_to_freq_unc=1e9*abs((mdl.Coefficients.Estimate(1)/mdl.Coefficients.Estimate(2)))*...
    sqrt((mdl.Coefficients.SE(1)/mdl.Coefficients.Estimate(1))^2+...
    (mdl.Coefficients.SE(2)/mdl.Coefficients.Estimate(2))^2);
new_to_wav_val=const.c/(new_to_freq_val*2);
new_to_wav_unc=new_to_freq_unc*const.c/((new_to_freq_val*2)^2);
fprintf('freq %.2f±%.2f MHz\n',new_to_freq_val*1e-6,new_to_freq_unc*1e-6)
fprintf('wavelength %.6f±%f nm \n',new_to_wav_val*1e9,new_to_wav_unc*1e9)
fprintf('error %e±%e nm \n',(new_to_wav_val-old_to_wav)*1e9,new_to_wav_unc*1e9)
fprintf('number of shots used %u\n',sum(data.mcp_tdc.num_ok))

toc

