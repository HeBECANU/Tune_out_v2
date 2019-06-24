function dum = amp_to_plot()
%%plot of the tune out value versus kick amplitude
%% setup directories
loop_config.dir = {
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190221_to_amp_3\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190221_to_amp_7\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190221_to_amp_10\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190214_to_hwp_168.5_nuller_reconfig_new_fiber_pdset_1.0v\'
    };


% Configs for plot
disp_config.num_bin = 4;
disp_config.colors_main=[[233,87,0];[33,188,44];[0,165,166]];
disp_config.plot_title='Amplitude dependence';
disp_config.font_name='cmr10';
disp_config.font_size_global=14;
disp_config.mdl_fun = @(b,x) b(1)+b(2).*x(:,1);
disp_config.beta0 = [1e14,1e5]; %Maybe check this
disp_config.opts = statset('nlinfit');
disp_config.fig_number = 3248;

% specify data
data = to_data(loop_config);
X=-data.drift.avg_coef(:,1);
Y=data.drift.to_val{1}.*1e-6;
Y_err=data.drift.to_val{2}.*1e-6;

% Plot the graph
plot_sexy(disp_config,X,Y,Y_err)
end
