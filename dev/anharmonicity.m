
anal_opts_osc_fit_cp=anal_opts.osc_fit;
anal_opts_osc_fit_cp.plot_fits=false;
anal_opts_osc_fit_cp.segment_borders = [1,16;
                                        17,34;
                                        35,55;
                                        56,78;
                                        79,100];
anal_opts_osc_fit_cp.num_seg = size(anal_opts_osc_fit_cp.segment_borders,1);
% data.mod_osc_fit = fit_trap_freq_mod(anal_opts_osc_fit_cp,data);

amps = cell(anal_opts_osc_fit_cp.num_seg,1);
freqs = cell(anal_opts_osc_fit_cp.num_seg,1);
% function anharmonic(data,anal_opts)
close all


osc_fit_master = data.mod_osc_fit;
rmse_old = ones(size(data.osc_fit.fit_rmse,1));
for segment = 1:anal_opts_osc_fit_cp.num_seg
    clear osc_fit
    osc_fit_seg.model_coefs = osc_fit_master.model_coefs(:,:,:,segment);
    osc_fit_seg.model = osc_fit_master.model(:,segment);
    osc_fit_seg.fit_rmse = osc_fit_master.fit_rmse(:,segment);
    
    osc_fit_seg.ok.did_fits=~cellfun(@(x) isequal(x,[]),osc_fit_seg.model);    
    mean_fit_rmse=nanmean(osc_fit_seg.fit_rmse(osc_fit_seg.ok.did_fits));
    median_fit_rmse=nanmedian(osc_fit_seg.fit_rmse(osc_fit_seg.ok.did_fits));
    std_fit_rmse=nanstd(osc_fit_seg.fit_rmse(osc_fit_seg.ok.did_fits));
    mean_fit_freq=nanmean(osc_fit_seg.model_coefs(osc_fit_seg.ok.did_fits,2,1));
    median_fit_freq=nanmedian(osc_fit_seg.model_coefs(osc_fit_seg.ok.did_fits,2,1));
    std_fit_freq=nanstd(osc_fit_seg.model_coefs(osc_fit_seg.ok.did_fits,2,1));
%     fprintf('fit error: median %f mean %f std %f\n',...
%        median_fit_rmse,mean_fit_rmse,std_fit_rmse)
    
    %label anything std_cut_fac*std away from the median as a bad fit
    std_cut_fac=1;
    freq_cut_fac = 3;
    mask=osc_fit_seg.ok.did_fits;
    osc_fit_seg.ok.rmse = osc_fit_seg.fit_rmse< median_fit_rmse+std_cut_fac*std_fit_rmse;
%     osc_fit_seg.ok.filter = osc_fit_seg.model_coefs(:,2)< mean_fit_freq + freq_cut_fac*std_fit_freq &...
%                             osc_fit_seg.model_coefs(:,2)> mean_fit_freq - freq_cut_fac*std_fit_freq; 
    %Removing annoying outliers by tricky insider knowledge
    osc_fit_seg.ok.filter = osc_fit_seg.model_coefs(:,2)< 60 &...
                            osc_fit_seg.model_coefs(:,2)> 40; 
                            
%     median_fit_rmse+std_cut_fac*std_fit_rmse
    calibration_shots = data.mcp_tdc.probe.calibration==1;
    measurement_shots = data.mcp_tdc.probe.calibration==0;
    osc_fit_seg.ok.cal=mask & osc_fit_seg.ok.rmse & osc_fit_seg.ok.filter & calibration_shots;
    osc_fit_seg.ok.meas=mask & osc_fit_seg.ok.rmse & osc_fit_seg.ok.filter & measurement_shots;
%     rmse_old = osc_fit_seg.fit_rmse;
    
    ok_amps_cal = osc_fit_seg.model_coefs(osc_fit_seg.ok.cal,1,1);
    ok_freqs_cal= osc_fit_seg.model_coefs(osc_fit_seg.ok.cal,2,1);
    ok_amps_meas = osc_fit_seg.model_coefs(osc_fit_seg.ok.meas,1,1);
    ok_freqs_meas= osc_fit_seg.model_coefs(osc_fit_seg.ok.meas,2,1);
    
    amps_cal{segment}  = abs(ok_amps_cal);
    amps_meas{segment} = abs(ok_amps_meas);
    freqs_cal{segment} = ok_freqs_cal; 
    freqs_meas{segment}= ok_freqs_meas;
    probe_freqs_meas{segment} = data.mcp_tdc.probe.freq.act.mean(osc_fit_seg.ok.meas);
    probe_freqs_cal{segment} = data.mcp_tdc.probe.freq.act.mean(osc_fit_seg.ok.cal);

end




amps_meas_all = (vertcat(amps_meas{:}));
freq_meas_all = (vertcat(freqs_meas{:}));
amps_cal_all = (vertcat(amps_cal{:}));
freq_cal_all = (vertcat(freqs_cal{:}));
probe_meas_all = vertcat(probe_freqs_meas{:});
probe_cal_all = vertcat(probe_freqs_cal{:});

num_slices = 15;
slice_edges = linspace(min(amps_meas_all),max(amps_meas_all),num_slices);
probe_bin_edges = linspace(min(probe_meas_all),max(probe_meas_all),probe_slices);

probe_trend = zeros(num_slices-1,4);
diff_trend = zeros(num_slices-1,4);
cal_trend = zeros(num_slices-1,4);

for i=1:num_slices-1
   c_bin_mask = amps_cal_all>=slice_edges(i) & amps_cal_all<=slice_edges(i+1);
   c_amps_bin = amps_cal_all(c_bin_mask);
   c_freqs_bin = freq_cal_all(c_bin_mask);
   N_c = sum(c_bin_mask);
   cal_trend(i,:) = [mean(c_amps_bin),mean(c_freqs_bin),std(c_amps_bin)/sqrt(N_c),std(c_freqs_bin)/sqrt(N_c)];
   
   p_bin_mask = amps_meas_all>=slice_edges(i) & amps_meas_all<=slice_edges(i+1);
   p_amps_bin = amps_meas_all(p_bin_mask);
   p_freqs_bin = freq_meas_all(p_bin_mask);
   N_p = sum(p_bin_mask);
   probe_trend(i,:) = [mean(p_amps_bin),mean(p_freqs_bin),std(p_amps_bin)/sqrt(N_p),std(p_freqs_bin)/sqrt(N_p)];
   diff_trend(i,:) = [mean(p_amps_bin),mean(p_freqs_bin) - cal_trend(i,2),...
       std(p_amps_bin)/sqrt(N_p),sqrt(std(p_freqs_bin)^2+cal_trend(i,4)^2)/sqrt(N_p)];
end


figure()

subplot(2,3,3)
num_freq_bins = 5;
probe_bin_edges = linspace(min(probe_meas_all),max(probe_meas_all),num_freq_bins);
amp_bin_edges = linspace(min(amps_meas_all),max(amps_meas_all),num_slices);
for i = 1:num_freq_bins-1
    freq_mask = probe_meas_all > probe_bin_edges(i) & probe_meas_all < probe_bin_edges(i+1) ;
    f_amp_bin = amps_meas_all(freq_mask);
    f_freq_bin = freq_meas_all(freq_mask);
    freq_slice = zeros(num_slices-1, 4);
    for j=1:num_slices-1
        amp_mask = f_amp_bin >= amp_bin_edges(j) & f_amp_bin < amp_bin_edges(j+1);
        fa_set = f_amp_bin(amp_mask);
        ff_set = f_freq_bin(amp_mask);
        N_f = numel(fa_set);
        freq_slice(j,:) = [mean(fa_set),mean(ff_set),std(fa_set)/sqrt(N_f), std(ff_set)/sqrt(N_f)];
    end
    c_f = viridis(num_freq_bins);
    e_f = errorbar(freq_slice(2:end,1),freq_slice(2:end,2)-cal_trend(i,2),freq_slice(2:end,4),freq_slice(2:end,4),freq_slice(2:end,3),freq_slice(2:end,3));
    hold on
    e_f.Color = c_f(i,:);
    xlabel('Osc amplitude')
    ylabel('Change in trap freq')
    title('Colour ~ probe freq')
%     colorbar()
end

%X = amp, Y = osc, bin by probe

%% figure drawing


colormap(viridis(1000))
subplot(2,3,1)
colordat = (probe_meas_all-min(probe_meas_all))/(max(probe_meas_all)-min(probe_meas_all));
scatter((abs(amps_meas_all)),(freq_meas_all),10,colordat,'filled')
title('Measurement fits')
xlabel('Peak amplitude')
ylabel('Fitted osc freq')


subplot(2,3,2)
colordat = (probe_cal_all-min(probe_cal_all))/(max(probe_cal_all)-min(probe_cal_all));
scatter((abs(amps_cal_all)),(freq_cal_all),10,colordat,'filled')
title('Calibration fits')
xlabel('Peak amplitude')
ylabel('Fitted osc freq')

% subplot(2,3,3)
% cf = viridis(probe_slices);
% colorinfo = (amps_meas_all-min(amps_meas_all))/(max(amps_meas_all)-min(amps_meas_all));
% ef = scatter(probe_meas_all,freq_meas_all,10,colorinfo,'filled');
% title('Freq change vs amplitude')
% xlabel('Probe freq')
% ylabel('Osc freq')

subplot(2,2,3)
et=errorbar(cal_trend(:,1),cal_trend(:,2),cal_trend(:,4),cal_trend(:,4),cal_trend(:,3),cal_trend(:,3),...
    '.','MarkerSize',20);
et.Color = cm(20,:);
hold on
pt=errorbar(probe_trend(:,1),probe_trend(:,2),probe_trend(:,4),probe_trend(:,4),probe_trend(:,3),probe_trend(:,3),...
    '.','MarkerSize',20);
pt.Color = cm(60,:);
legend('Probe off','Probe on')
title('Amplitude dependence')
xlabel('Peak amplitude')
ylabel('Fitted freq')
xlim([1,6])
ylim([47.5,48.6])

subplot(2,2,4)
ft=errorbar(diff_trend(:,1),diff_trend(:,2),diff_trend(:,4),diff_trend(:,4),diff_trend(:,3),diff_trend(:,3),...
    '.','MarkerSize',20);
hold on
cg = gray(20);
df=plot([1,6],[0,0]);
df.Color = cg(10,:);
xlim([1,6])
ylim([-0.1,0.5])
title('Fit freq difference')
xlabel('Oscillation amplitude')
ylabel('\omega_{test}-\omega_{cal} [Hz]')

suptitle('Comparison of anharmonic & probe effects')




