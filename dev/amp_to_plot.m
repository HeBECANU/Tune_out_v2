%%plot of the tune out value versus kick amplitude
clear all
%setup directories you wish to loop over
loop_config.dir = {
    'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20190221_to_amp_10\',
    'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20190221_to_amp_7\',
    'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20190221_to_amp_10\'
    };

data = to_data(loop_config);


weights = 1./(data.drift.to_val{2}.*1e6).^2;

sfigure(2909);
corr_plot(data.drift.avg_coef(:,1),data.drift.to_val{1},weights)
xlabel('Amp (Probe)')
ylabel(' Residual (MHz)')

sfigure(2910);
corr_plot(data.drift.avg_coef_cal(:,1),data.drift.to_val{1},weights)
xlabel('Amp (Cal)')
ylabel(' Residual (MHz)')