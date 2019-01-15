%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2018-08-14
% Author: J Ross
% A script to analyse the data we collected from a CW atom laser whose output flux is modulated by application of a probe beam near the tuneout wavelength.
% The script needs to:
%   * Imports all the data files in txy format
%   * Imports the log file, which indexes all the data as a) data run,
%       including wavelength, or b) calibration run to check the control
%       signal
%   * extracts the response at each set point, and
%   * Produces a plot of the measured response vs probe
%       wavelength, identifying the tuneout wavelength and giving a
%       stastistical uncertainty.

%Known BUGS/ possible Improvements--------------------------------------------------------
%2018-08-16
%histograming is done slightly incorrectly as the size of the bins is not excatly what you end up with
%fixed fft_tx
%made data structure for freq indexes

%Fixed !
%2018-08-16
%putting all the files in the apropriate main directory
%user var section seperate from code
%new import code used (faster)
%fixed histgraming error that was off by 1 bin
%fixed steady progression of wavelength assumption

%close all
%clear all
%%
% BEGIN USER VAR-------------------------------------------------

import_opts.dir='\\amplpc29\Users\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20180813_CW_AL_tuneout_scan\';
import_opts.file_name='d';
import_opts.force_reimport=0;
import_opts.force_forc=0;
import_opts.dld_xy_rot=0.61;
%Should probably try optimizing these
import_opts.xlim=[-20e-3, 20e-3];     %tight XY lims to eliminate hot spot from destroying pulse widths
import_opts.ylim=[-20e-3, 20e-3];
import_opts.tlim=[.4,1.4];
import_opts.txylim=[import_opts.tlim;import_opts.xlim;import_opts.ylim];

anal_opts.hist_tbin = 0.00011;
anal_opts.mod_freq=427;

% END USER VAR-----------------------------------------------------------


addpath('Colormaps') 
import_opts.shot_num=find_data_files(import_opts);
data=import_data(import_opts);

%now set up to histogram the time data
%general problem with histograming is when the time window is not an exact integer number of desired bin size wide
%the least bad solution is to round up the number of bins and have a half full bin at the end
num_shots = length(import_opts.shot_num);
num_bins = ceil(diff(import_opts.tlim)/anal_opts.hist_tbin); %round up 
hist_bin_edges = linspace(import_opts.tlim(1),import_opts.tlim(1)+num_bins*anal_opts.hist_tbin,num_bins);
%previously this assumed that the data was constant 
%hist_bin_centers = hist_bin_edges(1:end-1) + 0.5*t_stepsize;
hist_bin_centers = 0.5*(hist_bin_edges(1:end-1) + hist_bin_edges(2:end));


%will not work with data that has calibration shots
log.dir = strcat(import_opts.dir,'log_LabviewMatlab.txt');
log.table = readtable(log_dir,'Delimiter',',');
log.iter_nums = log_table.Var7;
log.setpoints = log_table.Var5*1e6;%convert to Hz
data.setpoints=log.setpoints;
log.psx_times = log_table.Var1;
set_freqs = unique(log.setpoints);


for ii=1:numel(set_freqs)
    freq=set_freqs(ii);
    %get the itteration numbers that had wavelength equaling freq
    matching_itt_num=log.iter_nums(freq==log.setpoints);
    data.freq_sorted.freq(ii)=freq;
    data.freq_sorted.indexes{ii}=matching_itt_num;
    %I think this is cleaner than copying data arround a lot, however if there is a lot of reporcessing going on the
    %line below can just write it out
    %data.freq_sorted.txy={data.txy{matching_itt_num}};
end
%usage example
%to get all the files that weer at freq =  data.freq_sorted.freq(5)
%{data.txy{data.freq_sorted.indexes{5}}}
%gives a cell array each of which contains matricies of a shots txy data
%vertcat(data.txy{data.freq_sorted.indexes{5}})
%gives all the data in a given wavelength

%%

plots=false;
for ii=1:size(data.txy,2)
    fprintf('%u\n',ii)
    hist_counts_raw=histcounts(data.txy{ii}(:,1),hist_bin_edges);
    hist_counts_smooth=gaussfilt(hist_bin_centers,hist_counts_raw,1e-4);
    fftout=fft_tx(hist_bin_centers,hist_counts_raw,10);
    [~,index]=min(abs(fftout(1,:)-anal_opts.mod_freq));
    %sample_freq=fftout(1,index)
    data.mod_amp_fft(ii)=fftout(2,index);
    if plots
        sfigure(1);
        clf
        subplot(1,2,1)
        plot(hist_bin_centers,hist_counts_smooth) 
        subplot(1,2,2)
        plot(fftout(1,:),abs(fftout(2,:)))
        xlim([400,450]) 
        pause(0.01)
    end
end

%%
sfigure(2);
subplot(1,2,1)
plot(data.setpoints,abs(data.mod_amp_fft),'x')
subplot(1,2,2)
plot(data.setpoints,angle(data.mod_amp_fft),'x')

%%
%find the mean angle at the start
indicies=data.setpoints<3.62855e14;
mean_angle=mean(angle(data.mod_amp_fft(indicies)))
theta = mean_angle; % to rotate 90 counterclockwise

%%


vector_rep=[real(data.mod_amp_fft);imag(data.mod_amp_fft)];
%dot(vector_rep(:,1),in_phase_unit)

R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
%create unit vectors
in_phase_unit = R*([1 0]');
out_phase_unit = R*([0 1]');

in_phase_comp=dot(repmat(in_phase_unit,[1,size(vector_rep,2)]),vector_rep,1);
out_phase_comp=dot(repmat(out_phase_unit,[1,size(vector_rep,2)]),vector_rep,1);

sfigure(3);
subplot(1,2,1)
plot(data.setpoints,in_phase_comp,'x')
subplot(1,2,2)
plot(data.setpoints,out_phase_comp,'x')



%%
cpx_rep=exp(theta*1j*pi).*data.mod_amp_fft;

sfigure(3);
subplot(1,2,1)
plot(data.setpoints,imag(cpx_rep),'x')
subplot(1,2,2)
plot(data.setpoints,real(cpx_rep),'x')

mean(angle(cpx_rep(indicies)))

%%
%the other way of doing this is more similar to what a lock in does
li_angle=-7.439045911031645e-01;
for ii=1:size(data.txy,2)
    fprintf('%u\n',ii)
    hist_counts_raw=histcounts(data.txy{ii}(:,1),hist_bin_edges);
    hist_counts_smooth=gaussfilt(hist_bin_centers,hist_counts_raw,1e-4);
    signal_in_phase=sum(hist_counts_smooth.*sin(2*pi*anal_opts.mod_freq*hist_bin_centers+li_angle))/(range(hist_bin_centers)*numel(hist_bin_centers));
    signal_out_phase=sum(hist_counts_smooth.*cos(2*pi*anal_opts.mod_freq*hist_bin_centers+li_angle))/(range(hist_bin_centers)*numel(hist_bin_centers));
    fftout=fft_tx(hist_bin_centers',hist_counts_smooth',10);
    data.mod_amp_li(ii,:)=[signal_out_phase,signal_in_phase];
end

%%
vector_rep=data.mod_amp_li';
%dot(vector_rep(:,1),in_phase_unit)
%theta = 7.439045911031645e-01;
theta=0;
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
%create unit vectors
in_phase_unit = R*([1 0]');
out_phase_unit = R*([0 1]');

in_phase_comp=dot(repmat(in_phase_unit,[1,size(vector_rep,2)]),vector_rep,1);
out_phase_comp=dot(repmat(out_phase_unit,[1,size(vector_rep,2)]),vector_rep,1);

sfigure(3);
subplot(1,2,1)
plot(data.setpoints,in_phase_comp,'x')
subplot(1,2,2)
plot(data.setpoints,out_phase_comp,'x')


%%

theta = 0.045; % to rotate 90 counterclockwise
li_data=data.mod_amp_li;
M=exp(-theta*1j).*(li_data(:,1)+li_data(:,2)*1j);
li_data = [real(M), imag(M)];

sfigure(4);
subplot(1,2,1)
plot(data.setpoints,li_data(:,1),'x')
subplot(1,2,2)
plot(data.setpoints,li_data(:,2),'x')


%%
%find the mean angle at the start
indicies=data.setpoints<3.6285e14;
li_angle=mean(angle(data.mod_amp_li(:,1)+data.mod_amp_li(:,2)*1j))
%%
% We know that the set point takes evenly spaced steps and is cyclic over
% 100 shots, so we can just say:
% wl_shots = cell(100,1); %Initialize cell for our data
% for i = 1:100
%     wl_shots{i} = []; %Create an empty array in each cell entry for easy use later
% end
% 
% for i=1:num_shots
%    idx = mod(i,100)+1;
%    wl_shots{idx} = [wl_shots{idx};data.txy{i}]; %Groups all shots of a particular probe wavelength into one compiled file.
% end
% 
% % Analyse bulk data
% com_sig=phasor_tranform(hist_bin_edges,data.txy,setpoints,427); %Produces complex FFT
% phasor_plot(com_sig,1) %Plots the outcome
% signed_signal = phasor_projection(com_sig,set_freqs,3.62855e8); %uses a phasor projection method to correct for the sign flip in response
% % This analysis removes outliers from each setpoint
% % The number of sd away from mean that are culled is set by thresh
% thresh = 1.5;
% trim_data = outlier_removal(signed_signal,thresh);
% 
% 
% % This combines all shots by wavlength setpoint, but the defective shots are
%  %detrimental here, so I suggest looking for another way...
% compressed_sig = phasor_tranform(hist_bin_edges,wl_shots,setpoints,427);
% signed_avg = phasor_projection(compressed_sig,set_freqs,3.62855e8);
% phasor_plot(compressed_sig,0)


%% Assuming good signals, let's set up a fit.
signal = trim_data;
fit_region = [3.6286,3.62875]*1e8;
mid_freq = median(fit_region);
fit_mask = find((signal(:,1)>fit_region(1))&(signal(:,1)<fit_region(2)));
fit_x = (signal(fit_mask,1)- mid_freq)/1e3; %Scales & centres region for better fit behavior
fit_y = signal(fit_mask,2);

%Tried a polyfit but couldn't get the root finder to work properly, this does fine
p_fun = @(p,x) p(4)+p(3)*x + p(2)*x.^2+p(1)*x.^3;
fit = fitnlm(fit_x,fit_y,p_fun,[0.5,-0.1,0.1,0]);
pc = fit.Coefficients.Estimate;
x = roots(pc);
x_real = x(imag(x)==0);


%Show me the money!
figure()
subplot(2,2,1)
scatter(signed_signal(:,1),signed_signal(:,2),'.');
title('Phase-corrected response')
xlabel('Frequency (MHz)')
ylabel('Modulation (kHz)')
subplot(2,2,2)
errorbar(trim_data(:,1),trim_data(:,2),trim_data(:,3))
title('trimmed data')
xlabel('Frequency (MHz)')
ylabel('Modulation (kHz)')
subplot(2,2,3)
plot(fit_x,fit_y,'.')
hold on
plot(fit_x,p_fun(fit.Coefficients.Estimate,fit_x));
xlabel(['\Delta=f - ',num2str(mid_freq),', GHz'])
ylabel('Modulation response')
title(['Zero crossing \Delta=',num2str(x_real),'GHz'])
subplot(2,2,4)
plot(p_fun(fit.Coefficients.Estimate,fit_x)-fit_y,'.')

%These errors are almost definitely wrong...
fit_errs = p_fun(fit.Coefficients.Estimate,fit_x)-fit_y;
SE = std(fit_errs);
MU = mean(fit_errs);

%Motivation: Find out where the curve crosses within 1sd of zero crossing,
%as this could be the ambiguity... Produces errors that seem too small to
%be believed.
x_pos_err = roots(pc-[0,0,0,SE]');
x_neg_err = roots(pc+[0,0,0,SE]');

format long
TOF = (2*mid_freq + x_real)*1e6 %TO frequency in blue
err_pos = 2*(x_real-x_pos_err(imag(x_pos_err)==0))*1e6
err_neg = 2*(x_real-x_neg_err(imag(x_neg_err)==0))*1e6
TOW = 299792458/TOF %TO wavelength in blue
W_err = -(299797458/TOF^2)*(abs(err_pos)+abs(err_neg))



%% Diagnosing spurious signals
    
%     f_test = 3.6286e8;
%     [~,t_idx] = min(abs(setpoints-f_test));
%     f_test = setpoints(t_idx);
%     test_shots = find(setpoints == f_test);
%     test_data = data.txy(test_shots);
%     figure()
%     for i=1:length(test_data)
%         ctr = ctr + 1;
%         temp_t= test_data{i}(:,1);
%         flux_t = histcounts(temp_t,hist_bin_edges);
%         [ft,c]= fft_tx(hist_bin_centers,flux_t');
%         f_c = 427;
%         [~,f_idx] = min(abs(ft(:,1)-427));
%         f_hw = 75;
%         [~,fw_idx] = min(abs(ft(:,1)-(f_c+f_hw)));
%         d_idx = fw_idx - f_idx;
%         f_window = [f_idx-d_idx:f_idx + d_idx];
%         signal = ft(f_window,2);     
%         plot(ft(f_window,1),signal)
%         hold on
%     end
%     title(['Response spectrum for f=',num2str(f_test),'Hz'])
%     xlabel('Frequency (Hz)')
%     ylabel('Response (kHz)')
    
    

    
   

