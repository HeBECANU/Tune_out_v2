clear all
%% Theory
to_freq_theory = 725736.480*1e9; %best theory guess in Hz
to_wav_theory = 299792458./to_freq_theory;
%% Exp parameters
aom_freq =0; %aom offset in Hz, as of 20190528 offset is calculated in main_trap_freq

%%
% find this .m file's path, this must be in the project root dir
this_folder = fileparts(which(mfilename));
% Add that folder plus all subfolders to the path.
addpath(genpath(this_folder));%add all subfolders to the path to find genpath_exclude
path_to_genpath=fileparts(which('genpath_exclude'));
path(pathdef) %clean up the path back to the default state to remove all the .git that were added
addpath(this_folder)
addpath(path_to_genpath)
addpath(genpath_exclude(fullfile(this_folder,'lib'),'\.')) %dont add hidden folders
addpath(genpath_exclude(fullfile(this_folder,'dev'),'\.'))
addpath(genpath_exclude(fullfile(this_folder,'bin'),'\.'))

%% polarisation data options
pol_opts.location = 'pre_right';%post, pre_cen, pre_left, pre_right
pol_opts.predict = 'fit';%'interp'; %obs (obsovation) fit (pertial fit) full_fit (fit with all parameters free)
%% Display options
opts = statset('MaxIter',1e4);
ci_size=1-erf(1/sqrt(2));%1-erf(zvalue/sqrt(2)) %confidence interval for cutting outliers
disp_config.font_name = 'Times New Roman';%'cmr10';
disp_config.font_size_global=14;
disp_config.opts=statset('nlinfit');
disp_config.colors_main = [[75,151,201];[193,114,66];[87,157,95]];
disp_config.bin_tol=0.01;
%% HWP
loop_config.dir = {
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190204_to_hwp_51_nuller_reconfig\'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190203_to_hwp_29_nuller_reconfig\'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190201_to_hwp_131_nuller_reconfig\'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190203_to_hwp_171_nuller_reconfig\'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190203_to_hwp_191_nuller_reconfig\'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190211_to_hwp_208_nuller_reconfig\'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\\20190211_to_hwp_194_nuller_reconfig\'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190210_to_hwp_187_nuller_reconfig\'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190210_to_hwp_181_nuller_reconfig\'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190210_to_hwp_177_nuller_reconfig_part_b\'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190210_to_hwp_177_nuller_reconfig_part_a\'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190210_to_hwp_171_nuller_reconfig\'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190210_to_hwp_165_nuller_reconfig\'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190210_to_hwp_160_nuller_reconfig\'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190209_to_hwp_155_nuller_reconfig\'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190209_to_hwp_145_nuller_reconfig\'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190209_to_hwp_140_nuller_reconfig\'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190208_to_hwp_120_nuller_reconfig_okish\'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190208_to_hwp_99_nuller_reconfig\'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190207_to_hwp_80_nuller_reconfig\'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190204_to_hwp_121_nuller_reconfig\'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190204_to_hwp_24_nuller_reconfig\'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190204_to_hwp_46_nuller_reconfig\'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190204_to_hwp_61_nuller_reconfig\'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190204_to_hwp_92_nuller_reconfig\'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190211_to_hwp_230_nuller_reconfig\'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190211_to_hwp_217_nuller_reconfig\'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190211_to_hwp_199_nuller_reconfig\'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190206_to_hwp_100_nuller_reconfig\'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190204_to_hwp_141_nuller_reconfig\'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190204_to_hwp_160_nuller_reconfig\'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190204_to_hwp_180_nuller_reconfig\'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190204_to_hwp_70_nuller_reconfig\'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190211_to_hwp_240_nuller_reconfig\'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190211_to_hwp_250_nuller_reconfig_tenma_setpoint\'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190201_to_hwp_111_nuller_reconfig\'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190203_to_hwp_151_nuller_reconfig\'
};
data.hwp = to_data(loop_config);
%% QWP
loop_config.dir = {
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190304_qwp_280_pure'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190302_qwp_283_pure_long'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190301_qwp_283_pure_run'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190301_qwp_246_2'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190227_qwp_270'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190227_qwp_286'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190227_qwp_310'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190227_qwp_286_no_analog'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190226_qwp_254'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190226_qwp_246'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190226_qwp_234'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190225_qwp_226'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190225_qwp_220'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190225_qwp_202'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190225_qwp_187'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190225_qwp_177'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190225_qwp_162'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190225_qwp_154'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190224_qwp_130'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190224_qwp_134'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190224_qwp_138'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190224_qwp_142'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190224_qwp_146'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190223_qwp_150'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190226_qwp_260'
};
data.qwp = to_data(loop_config);
%% MIX
loop_config.dir = {
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190305_hwp_340_qwp_270'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190306_hwp_350_qwp_274'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190306_overnight_hwp_354_qwp_268'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190307_hwp_10_qwp_280'
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190307_hwp_06_qwp_290'
};
data.mix = load_pocessed_to_data(loop_config);
%% Generate data vectors
to_val = [data.hwp.drift.to_val{1};data.qwp.drift.to_val{1};data.mix.drift.to_val{1}];%all our single scan tune-out values
to_unc = [data.hwp.drift.to_val{2};data.qwp.drift.to_val{2};data.mix.drift.to_val{2}];%all corresponding uncertainties
wlin=1./(to_unc.^2);

to_val_run = [data.hwp.main.lin_fit{1};data.qwp.main.lin_fit{1};data.mix.main.lin_fit{1}];%the to val for a particular run

V = zeros(numel(to_val),1);%fourth stokes parameter
d_p = zeros(numel(to_val),1);%power diffrence between min and max transmition
theta = zeros(numel(to_val),1);%min angle

[V_run,d_p_run,theta_run] = pol_data_import(pol_opts);

shot_idx = 1;
scan_num_vec = [data.hwp.main.scan_num;data.qwp.main.scan_num;data.mix.main.scan_num];

for loop_idx=1:size(V_run,1)
    scan_num = scan_num_vec(loop_idx);
    V(shot_idx:(shot_idx+scan_num-1),1) = ones(scan_num,1).*V_run(loop_idx);
    d_p(shot_idx:(shot_idx+scan_num-1),1) = ones(scan_num,1).*d_p_run(loop_idx);
    theta(shot_idx:(shot_idx+scan_num-1),1) = ones(scan_num,1).*theta_run(loop_idx);
    shot_idx = shot_idx+scan_num;
end

%% Fit the full model to all of our data
Q_fun = @(d_p,theta,phi) d_p.*cos(2.*(theta+phi.*pi));%function to calculate second stokes parameter
full_mdl = @(b,x) b(1) + b(2).*x(:,2) + b(3).*(1+Q_fun(x(:,1),x(:,3),b(4)));%full model for how the tune out behaves
vars = [d_p,V,theta];
to_val_fit = to_val;
wlin_fit = wlin./sum(wlin);
beta0 = [7.257355*1e14,6*1e9,3*1e9,0.8];%0.2 or 0.8
fit_mdl = fitnlm(vars,to_val_fit,full_mdl,beta0,...
    'Options',opts,'Weights',wlin_fit,'CoefficientNames' ,{'tune_out','vector','tensor','phase'});
fit_vals = fit_mdl.Coefficients{:,1};

%second fit with phase fixed
Q = Q_fun(d_p,theta,fit_vals(4));%second stokes parameter
Q_run = Q_fun(d_p_run,theta_run,fit_vals(4));%second stokes parameter for each run
full_mdl = @(b,x) b(1) + b(2).*x(:,2) + b(3).*(1+x(:,1));%full model for how the tune out behaves
vars = [Q,V];
beta0 = fit_vals(1:3);
fit_mdl = fitnlm(vars,to_val_fit,full_mdl,beta0,...
    'Options',opts,'Weights',wlin_fit,'CoefficientNames' ,{'tune_out','vector','tensor'});
fit_vals = fit_mdl.Coefficients{:,1};
fit_uncs = fit_mdl.Coefficients{:,2};

%% Full 2D plot in stokes space
% sfigure(1234);
% clf
% [R, T] = meshgrid(linspace(0,1),linspace(0,2*pi));
% TO = fit_vals(1)+fit_vals(2).*R.*cos(T)+fit_vals(3).*R.*sin(T);
% h = surface(R.*cos(T),R.*sin(T),TO./1e6-fit_vals(1)./1e6);
% colormap(viridis)
% set(h, 'FaceAlpha', 0.7)
% shading interp
% hold on
% scatter3(V,Q,to_val./1e6-fit_vals(1)./1e6,'MarkerFaceColor',[0 .75 .75])
% grid on
% title('Full stokes space representation')
% xlabel('Fourth Stokes Parameter, V')
% ylabel('Second Stokes Parameter, Q')
% zlabel( sprintf('Tune out value-%.1f±%.1f (MHz)',fit_vals(1)./1e6,fit_uncs(1)./1e6) )
% set(gcf,'color','w')
% set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1600, 900])

%% Plots of tune-out dependance on the individual stokes parameters
modelfun = @(b,x) b(1) + b(2).*x(:,1);

vec_corr_to = to_val-V.*fit_vals(2);
beta0 = [fit_vals(1),fit_vals(3)];
fit_mdl_t = fitnlm(Q,vec_corr_to./1e6,modelfun,beta0,...
    'Options',opts,'Weights',wlin,'CoefficientNames' ,{'tune_out','tensor'});

tens_corr_to = to_val-Q.*fit_vals(3);
beta0 = [fit_vals(1),fit_vals(2)];
fit_mdl_v = fitnlm(V,tens_corr_to./1e6,modelfun,beta0,...
    'Options',opts,'Weights',wlin,'CoefficientNames' ,{'tune_out','vector'});

disp_config.plot_title = '';
disp_config.x_label = 'Fourth Stokes Parameter, V';
disp_config.fig_number=11235;
disp_config.plot_offset.val=fit_vals(1)./1e6;
disp_config.plot_offset.unc=fit_uncs(1)./1e6;
plot_sexy(disp_config,V,tens_corr_to./1e6,wlin,fit_mdl_v)

disp_config.plot_title = '';
disp_config.x_label = 'Second Stokes Parameter, Q';
disp_config.fig_number=11236;
disp_config.plot_offset.val=fit_vals(1)./1e6;
disp_config.plot_offset.unc=fit_uncs(1)./1e6;
plot_sexy(disp_config,Q,vec_corr_to./1e6,wlin,fit_mdl_t)

%% Residuals
to_res_run = to_val_run-predict(fit_mdl,[Q_run,V_run]);
to_res = to_val-predict(fit_mdl,[Q,V]);
sfigure(33099);
plot(to_res_run./1e6)
xlabel('index')
ylabel('residual (Mhz)')
sfigure(4445005);
histfit(to_res./1e6,50)
xlabel('residual (MHz)')
ylabel('count')
%% Write out the results

to_freq_val=fit_vals(1);
to_freq_unc=fit_uncs(1);
to_wav_val=299792458/(to_freq_val);
to_wav_unc=2*to_freq_unc*299792458/(to_freq_val^2);
old_to_wav=413.0938e-9;

fprintf('\n====TO full fit results==========\n')
fprintf('TO freq             %.1f±(%.0f) MHz\n',...
    to_freq_val*1e-6,to_freq_unc*1e-6)
fprintf('TO wavelength       %.6f±%f nm \n',to_wav_val*1e9,to_wav_unc*1e9)
fprintf('diff from TOV1               %e±%e nm \n',(to_wav_val-old_to_wav)*1e9,to_wav_unc*1e9)
fprintf('diff from Theory             %e±%e nm \n',(to_wav_val-to_wav_theory)*1e9,to_wav_unc*1e9)