function disp_to_results(data,anal_opts)

global const

to_res_fit_all=data.to_fit_all;
to_fit_all_trimed_val={to_res_fit_all.fit_trimmed.to_freq(:).val};
to_fit_all_trimed_unc=to_res_fit_all.fit_trimmed.to_unc_boot;
to_fit_unc_unc_boot_lin=to_res_fit_all.fit_trimmed.boot{1}.results.se_se_fun_whole;
to_fit_unc_unc_boot_quad=to_res_fit_all.fit_trimmed.boot{2}.results.se_se_fun_whole;

to_fit_seg_val=data.to_fit_seg.fit_trimmed.to_freq.sewm.val;
to_fit_seg_unc=data.to_fit_seg.fit_trimmed.to_freq.sewm.unc;

%inverse scaled gradient to give the single shot uncert (with scaling factor to include calibration)
tot_num_shots=to_res_fit_all.num_shots+data.cal.num_shots;
single_shot_uncert=data.to_fit_all.fit_trimmed.single_shot_unc{1}...
    *sqrt(tot_num_shots/to_res_fit_all.num_shots);
fprintf('\n====TO fit results==========\n')
fprintf('dir =%s\n',anal_opts.tdc_import.dir)
%calculate some statistics and convert the model parameter into zero crossing and error therin
old_to_wav=413.0938e-9;

%to_res_fit_all.fit_trimmed.to_unc_fit
to_wav_val_lin=const.c/(to_fit_all_trimed_val{1});
to_freq_val_lin=to_fit_all_trimed_val{1};
to_freq_unc_lin=to_fit_all_trimed_unc{1};
to_wav_unc_lin=2*to_fit_all_trimed_unc{1}*const.c/((to_fit_all_trimed_val{1}*2)^2);
to_wav_val_quad=const.c/(to_fit_all_trimed_val{2});
to_freq_val_quad=to_fit_all_trimed_val{2};
to_freq_unc_quad=to_fit_all_trimed_unc{2};
to_wav_unc_quad=2*to_fit_all_trimed_unc{2}*const.c/((to_fit_all_trimed_val{2}*2)^2);
time_run_start=data.mcp_tdc.time_create_write(1,2)-anal_opts.trig_dld-anal_opts.dld_aquire;
time_run_stop=data.mcp_tdc.time_create_write(end,2)-anal_opts.trig_dld-anal_opts.dld_aquire;
time_duration=data.mcp_tdc.time_create_write(end,2)-data.mcp_tdc.time_create_write(1,2);
fprintf('run start time               %.1f        (posix)\n',time_run_start)
fprintf('                             %s (ISO)\n',datestr(datetime(time_run_start,'ConvertFrom','posix'),'yyyy-mm-ddTHH:MM:SS'))
fprintf('run stop time                %.1f        (posix)\n',time_run_stop)
fprintf('                             %s (ISO)\n',datestr(datetime(time_run_stop,'ConvertFrom','posix'),'yyyy-mm-ddTHH:MM:SS'))
fprintf('duration                     %.1f (s),%.1f (h)\n',...
    time_duration,time_duration/(60*60))

fprintf('median damping time %.2f\n',median(1./data.osc_fit.model_coefs(data.osc_fit.ok.rmse,7,1)))

fprintf('TO freq (Linear,all)         %.1f±(%.0f±%.0f) MHz\n',...
    to_freq_val_lin*1e-6,to_freq_unc_lin*1e-6,to_fit_unc_unc_boot_lin*1e-6*2)
fprintf('TO freq (Quadratic,all)      %.1f±(%.0f±%.0f) MHz\n',...
    to_freq_val_quad*1e-6,to_freq_unc_quad*1e-6,to_fit_unc_unc_boot_quad*1e-6*2)
fprintf('TO freq (Linear,seg)         %.1f±(%.0f) MHz\n',...
    to_fit_seg_val*1e-6,to_fit_seg_unc*1e-6)

fprintf('diff between Lin and Quad    %e±%e MHz \n',(to_freq_val_lin-to_freq_val_quad)*1e-6,sqrt(to_freq_unc_lin^2+to_freq_unc_quad^2)*1e-6)
fprintf('diff between all and seg     %e±%e MHz \n',(to_freq_val_lin-to_fit_seg_val)*1e-6,sqrt(to_freq_unc_lin^2+to_fit_seg_unc^2)*1e-6)

fprintf('TO wavelength (Linear)       %.6f±%f nm \n',to_wav_val_lin*1e9,to_wav_unc_lin*1e9)
fprintf('TO wavelength (Quadratic)    %.6f±%f nm \n',to_wav_val_quad*1e9,to_wav_unc_quad*1e9)
fprintf('diff between Lin and Quad    %e±%e nm \n',(to_wav_val_lin-to_wav_val_quad)*1e9,sqrt(to_wav_unc_lin^2+to_wav_unc_quad^2)*1e9)
fprintf('diff from TOV1               %e±%e nm \n',(to_wav_val_lin-old_to_wav)*1e9,to_wav_unc_lin*1e9)
%more logic needs to be included here
fprintf('files with enough number     %u\n',sum(data.mcp_tdc.num_ok'))
fprintf('number of probe files        %u \n',to_res_fit_all.num_shots)
fprintf('number of calibration files  %u \n',data.cal.num_shots)
fprintf('total used                   %u \n',tot_num_shots)
fprintf('shot uncert scaling @1SD %.1f MHz, %.2f fm /sqrt(shots)\n',single_shot_uncert*1e-6,...
    single_shot_uncert*const.c/((to_fit_all_trimed_val{1}*2)^2)*10^15)
%predicted uncert using this /sqrt(n), unless derived differently this is pointless
%fprintf('predicted stat. uncert %.1f MHz, %.2f fm\n',single_shot_uncert/sqrt(tot_num_shots)*1e-6,...
%    single_shot_uncert/sqrt(tot_num_shots)*const.c/((to_fit_trimed_val*2)^2)*10^15)
end