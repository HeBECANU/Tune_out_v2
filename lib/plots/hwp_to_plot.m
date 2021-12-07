function dum = hwp_to_plot()
%complete plot of hwp
%Script that scrapes the analysed data from dirs (currently messy but works)
%setup directories you wish to loop over
% BEGIN USER VAR-------------------------------------------------
root_dir = 'E:\Data\Tuneout\to_main_data';
all_folders = dir(root_dir);
all_folders = all_folders(3:end);
use_folders = all_folders(2:38);

loop_config.dir = {use_folders.name};
loop_config.dir = cellfun(@(x) fullfile(root_dir,x),loop_config.dir,'uni',0);
% END US
% END USER VAR-----------------------------------------------------------

%add all subfolders to the path
% this_folder = fileparts(which(mfilename));
% folder=strsplit(this_folder,filesep); %spllit path into folders
% folder=strjoin(folder(1:end-2),filesep); %go up two
% Add that folder plus all subfolders to the path.
% addpath(genpath(folder));


data = load_pocessed_to_data(loop_config);
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
    fin = strfind(current_dir,'_p');
    to_pol(loop_idx) = str2num(current_dir(strt+4:fin-1));
    to_pol_drift(shot_idx:(shot_idx+data.main.scan_num(loop_idx)-1),1) = ones(data.main.scan_num(loop_idx),1).*to_pol(loop_idx);
    shot_idx = shot_idx+data.main.scan_num(loop_idx);
end
% %
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
% <<<<<<< HEAD
disp_config.bin_tol=0.01;

modelfun = @(b,x) b(1).*(cos(x(:,1).*pi./180+b(2).*2*pi).^2)+b(3);
opts = statset('MaxIter',1e4);
weights = 1./to_vals_error.^2;
fit_mdl = fitnlm(to_pol_drift',vec_corr_to',modelfun,beta0,'Options',opts,'Weight',weights);%'ErrorModel','combined'
plot_sexy(disp_config,to_pol_drift,vec_corr_to,weights,fit_mdl)
end

