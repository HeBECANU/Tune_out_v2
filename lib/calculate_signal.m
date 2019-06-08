function data_signal=calculate_signal(anal_opts_fit_to,data)
%calculate the probe beam signal

%%input masking
is_cal=col_vec(data.mcp_tdc.probe.calibration); %because its used a lot make a temp var for calibration logic vector
is_cal(isnan(is_cal))=1;    
probe_dat_mask=col_vec(data.osc_fit.ok.all) & ~is_cal &  ~isnan(col_vec(data.wm_log.proc.probe.freq.act.mean))...
    & ~isnan(col_vec(data.osc_fit.trap_freq_recons));

%%selecting the data
trap_freq=col_vec(data.osc_fit.trap_freq_recons(probe_dat_mask));
trap_freq_unc=col_vec(data.osc_fit.trap_freq_recons_unc(probe_dat_mask));
cal_trap_freq=col_vec(data.cal.freq_drift_model(data.mcp_tdc.time_create_write(probe_dat_mask,2)));
cal_trap_freq_unc=col_vec(data.cal.unc);
if ~isequal(size(cal_trap_freq_unc),[1,1])
    error('not ready for cal unc vec')
end

%% probe beam frequency calculation
square_trap_freq= (trap_freq).^2-(cal_trap_freq).^2;
square_trap_freq_unc=sqrt(2).*sqrt((trap_freq_unc.*trap_freq).^2+(cal_trap_freq_unc.*cal_trap_freq).^2);

%% output the data
data_signal=[];
data_signal.square_probe_trap_freq.val=nan*is_cal;
data_signal.square_probe_trap_freq.unc=nan*is_cal;
data_signal.square_probe_trap_freq.val(probe_dat_mask)=square_trap_freq;
data_signal.square_probe_trap_freq.unc(probe_dat_mask)=square_trap_freq_unc;
data_signal.dat_mask=col_vec(probe_dat_mask);



end