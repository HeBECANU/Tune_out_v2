function dum = hwp_to_plot()
%complete plot of hwp
%Script that scrapes the analysed data from dirs (currently messy but works)
%setup directories you wish to loop over
% BEGIN USER VAR-------------------------------------------------
loop_config.dir = {
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190204_to_hwp_51_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190203_to_hwp_29_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190201_to_hwp_131_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190203_to_hwp_171_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190203_to_hwp_191_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190211_to_hwp_208_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\\20190211_to_hwp_194_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190210_to_hwp_187_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190210_to_hwp_181_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190210_to_hwp_177_nuller_reconfig_part_b\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190210_to_hwp_177_nuller_reconfig_part_a\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190210_to_hwp_171_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190210_to_hwp_165_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190210_to_hwp_160_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190209_to_hwp_155_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190209_to_hwp_145_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190209_to_hwp_140_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190208_to_hwp_120_nuller_reconfig_okish\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190208_to_hwp_99_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190207_to_hwp_80_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190204_to_hwp_121_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190204_to_hwp_24_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190204_to_hwp_46_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190204_to_hwp_61_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190204_to_hwp_92_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190211_to_hwp_230_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190211_to_hwp_217_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190211_to_hwp_199_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190206_to_hwp_100_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190204_to_hwp_141_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190204_to_hwp_160_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190204_to_hwp_180_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190204_to_hwp_70_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190211_to_hwp_240_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190211_to_hwp_250_nuller_reconfig_tenma_setpoint\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190201_to_hwp_111_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190203_to_hwp_151_nuller_reconfig\'
    };
% END USER VAR-----------------------------------------------------------

%add all subfolders to the path
this_folder = fileparts(which(mfilename));
folder=strsplit(this_folder,filesep); %spllit path into folders
folder=strjoin(folder(1:end-2),filesep); %go up two
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));


data = to_data(loop_config);
drift_data_compiled = data.drift;
main_data_compiled = data.main;
TO_st_pt = 7.257355*1e14;
selected_dirs = 1:numel(loop_config.dir);
to_pol = zeros(numel(loop_config.dir),1);
to_pol_drift = zeros(numel(data.drift.to_time),1);
shot_idx = 1;

for loop_idx=selected_dirs
    current_dir = loop_config.dir{loop_idx};
    strt = strfind(current_dir,'hwp_');
    fin = strfind(current_dir,'_n');
    to_pol(loop_idx) = str2num(current_dir(strt+4:fin-1));
    to_pol_drift(shot_idx:(shot_idx+data.main.scan_num(loop_idx)-1),1) = ones(data.main.scan_num(loop_idx),1).*to_pol(loop_idx);
    shot_idx = shot_idx+data.main.scan_num(loop_idx);
end
%%
vec_corr_to = data.drift.to_val{1}./1e6;
to_vals_error = data.drift.to_val{2}./1e6;

%fit sin waves to the two sets of data

beta0 = [1e3,0.5,nanmean(vec_corr_to)]; %intial guesses
disp_config.colors_main = [[233,87,0];[33,188,44];[0,165,166]];
disp_config.plot_title = 'Tune-out Dependence on Input Polarization Angle';
disp_config.x_label='Input Polarization angle (degrees)';
disp_config.font_name = 'cmr10';
disp_config.font_size_global=14;
disp_config.opts=statset('nlinfit');
disp_config.fig_number=3400;
<<<<<<< HEAD
disp_config.bin_tol=0.01;

modelfun = @(b,x) b(1).*(cos(x(:,1).*pi./180+b(2).*2*pi).^2)+b(3);
opts = statset('MaxIter',1e4);
weights = 1./to_vals_error.^2;
fit_mdl = fitnlm(to_pol_drift',vec_corr_to',modelfun,beta0,'Options',opts,'Weight',weights);%'ErrorModel','combined'
plot_sexy(disp_config,to_pol_drift,vec_corr_to,weights,fit_mdl)
end

