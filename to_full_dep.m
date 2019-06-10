clear all
%% Theory
to_freq_theory = 725736480*1e6; %best theory guess in Hz
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

%%
hebec_constants %call the constants function that makes some globals
%% polarisation data options
pol_opts.location = 'post';%post, pre_cen, pre_left, pre_right
pol_opts.predict = 'full_fit';%'full_fit';%'interp'; %obs (obsovation) fit (pertial fit) full_fit (fit with all parameters free)
%% Display options
opts = statset('MaxIter',1e4);
ci_size=1-erf(1/sqrt(2));%1-erf(zvalue/sqrt(2)) %confidence interval for cutting outliers
disp_config.font_name = 'Times New Roman';%'cmr10';
disp_config.font_size_global=14;
disp_config.opts=statset('nlinfit');
disp_config.colors_main = [[75,151,201];[193,114,66];[87,157,95]];
disp_config.bin_tol=0.01;
%% HWP
% loop_config.dir = {
%     '..\scratch_data\20190227_qwp_270',
%     '..\scratch_data\20190227_qwp_286',
%     '..\scratch_data\20190227_qwp_310',
%     };
%root_data_dir='..\scratch_data';
root_data_dir='G:\good_data';
files = dir(root_data_dir);
files=files(3:end);
% Get a logical vector that tells which is a directory.
dir_mask = [files.isdir];
folders=files(dir_mask);
folders=arrayfun(@(x) fullfile(root_data_dir,x.name),folders,'UniformOutput' ,false);
folders
loop_config.dir=folders;
%%

data = load_pocessed_to_data(loop_config);

%%
%save('20190608_imported_data_for_to_full_dep.mat')
load('20190608_imported_data_for_to_full_dep.mat')


%% Generate data vectors
%to_val = [data.hwp.drift.to_val{1};data.qwp.drift.to_val{1};data.mix.drift.to_val{1}];%all our single scan tune-out values
%to_unc = [data.hwp.drift.to_val{2};data.qwp.drift.to_val{2};data.mix.drift.to_val{2}];%all corresponding uncertainties
to_val_for_polz=data.drift.to.val;
to_unc=data.drift.to.unc;
wlin=1./(to_unc.^2);


pol_opts.hwp=data.drift.wp.hwp;
pol_opts.qwp=data.drift.wp.qwp;
pol_model=pol_data_query(pol_opts);

polz_theta=pol_model.theta.val;
polz_v=pol_model.v.val;
polz_cont=pol_model.cont.val;

%% Fit the full model to all of our data
Q_fun = @(d_p,theta,phi) d_p.*cos(2.*(theta+phi.*pi));%function to calculate second stokes parameter, in a rotated frame
full_mdl = @(b,x) b(1) + b(2).*x(:,2) + b(3).*(1+Q_fun(x(:,1),x(:,3),b(4)));%full model for how the tune out behaves
vars = [polz_cont,polz_v,polz_theta];
to_val_fit = to_val_for_polz;
wlin_fit = wlin./sum(wlin);
to_fit_val_offset=wmean(to_val_fit,wlin_fit);

 gf_opt=[];
gf_opt.domain=[[-1,1]*1e5;...   %offset
               [-1,1]*1e4;...  %vector
               [-1,1]*1e4;...  %tensor
               [-1,1]*4*pi;...  %phase
               ];        
gf_opt.start=[0, 1e3, 1*1e3, (rand(1)-0.5)*2*pi];
gf_opt.rmse_thresh=2e-3;
gf_opt.plot=false;
gf_opt.level=2;
gf_out=global_fit(vars,(to_val_fit-to_fit_val_offset)*1e-6,full_mdl,gf_opt);


%
opts = statset('MaxIter',1e4);
beta0 = gf_out.params%0.2 or 0.8
fit_mdl_full = fitnlm(vars,(to_val_fit-to_fit_val_offset)*1e-6,full_mdl,beta0,...
    'Options',opts,'Weights',wlin_fit,'CoefficientNames' ,{'tune_out','vector','tensor','phase'})
fit_vals_full = fit_mdl_full.Coefficients.Estimate;
%
%second fit with phase fixed
polz_q = Q_fun(polz_cont,polz_theta,fit_vals_full(4));%second stokes parameter, relative to the vector that is otrhogonal to the B feild and beam axis
% TO CHECK
%Q_run = Q_fun(polz_cont,polz_theta,fit_vals(4));%second stokes parameter for each run
full_mdl = @(b,x) b(1) + b(2).*x(:,2) + b(3).*(1+x(:,1));%full model for how the tune out behaves
vars = [polz_q,polz_v];
beta0 = fit_vals_full(1:3);
opts = statset('MaxIter',1e4);
fit_mdl = fitnlm(vars,(to_val_fit-to_fit_val_offset)*1e-6,full_mdl,beta0,...
    'Options',opts,'Weights',wlin_fit,'CoefficientNames' ,{'tune_out','vector','tensor'})
fit_coefs_fixed_phase.val = fit_mdl.Coefficients.Estimate;
fit_coefs_fixed_phase.unc = fit_mdl.Coefficients.SE;

to_scalar_minus_half_tensor=[];
to_scalar_minus_half_tensor.val=(fit_coefs_fixed_phase.val(1)-0*fit_coefs_fixed_phase.val(3))*1e6+to_fit_val_offset;
to_scalar_minus_half_tensor.unc=fit_coefs_fixed_phase.unc(1)*1e6;
%to_scalar_minus_half_tensor.val=(fit_coefs_fixed_phase.val(1)+1*fit_coefs_fixed_phase.val(3))*1e6+to_fit_val_offset;
%to_scalar_minus_half_tensor.unc=sqrt(fit_coefs_fixed_phase.unc(1)^2+(0.5*fit_coefs_fixed_phase.unc(3))^2)*1e6;
fprintf('%.1f\n',(to_scalar_minus_half_tensor.val*1e-6-725735000))


% 2   2000
% 0.5 1000
% 2/3 700
% Write out the results
%%
full_data_set=[polz_cont,polz_v,polz_theta,to_val_for_polz,to_unc];
full_data_set=num2cell(full_data_set,2);
% check that function gives the same answer as the above procedure
fprintf('%.1f\n',(two_stage_two_dim_to_fit(full_data_set)*1e-6-725735000))

%%


boot=bootstrap_se(@two_stage_two_dim_to_fit,full_data_set,...
    'plots',true,...
    'replace',true,...
    'samp_frac_lims',[0.01,1],...%[0.005,0.9]
    'num_samp_frac',100,... %20
    'num_samp_rep',100,... %1e2
    'plot_fig_name','TO fit bootstrap',...
    'save_multi_out',0,...
    'verbose',3);


to_scalar_minus_half_tensor.unc_boot=boot.results.se_fun_whole*1e6;
to_scalar_minus_half_tensor.unc_unc_boot=boot.results.se_se_fun_whole*1e6;

%%

to_wav_val=f2wl(to_scalar_minus_half_tensor.val);
to_wav_unc=2*to_scalar_minus_half_tensor.unc*const.c/(to_scalar_minus_half_tensor.val^2);
old_to_wav=413.0938e-9;

fprintf('\n====TO full fit results==========\n')
fprintf('TO freq             %.1f±(%.0f(fit),%.0f±%.0f(boot)) MHz\n',...
    to_scalar_minus_half_tensor.val*1e-6,...
    to_scalar_minus_half_tensor.unc*1e-6,...
    to_scalar_minus_half_tensor.unc_boot*1e-6,...
    to_scalar_minus_half_tensor.unc_unc_boot*1e-6)
fprintf('diff from TOV1      %.1f±%.1f MHz \n',(to_scalar_minus_half_tensor.val-f2wl(old_to_wav))*1e-6,to_scalar_minus_half_tensor.unc*1e-6)
fprintf('diff from Theory    %e±%e MHz \n',(to_scalar_minus_half_tensor.val-to_freq_theory)*1e-6,to_scalar_minus_half_tensor.unc*1e-6)
fprintf('TO wavelength       %.6f±%f nm \n',to_wav_val*1e9,to_wav_unc*1e9)

%% Full 2D plot in stokes space
stfig('2d TO fit');
clf
[surf_polz_q, surf_polz_v] = meshgrid(linspace(-1,1),linspace(-1,1));
%full_mdl = @(b,x) b(1) + b(2).*x(:,2) + b(3).*(1+Q_fun(x(:,1),x(:,3),b(4)));%full model for how the tune out behaves
%TO = fit_vals(1)+fit_vals(2).*R.*cos(T)+fit_vals(3).*R.*sin(T);
%h = surface(R.*cos(T),R.*sin(T),TO./1e6-fit_vals(1)./1e6);
Q_fun = @(d_p,theta,phi) d_p.*cos(2.*(theta+phi.*pi));%function to calculate second stokes parameter
surf_mdl = @(b,x) (b(1) + b(2).*x(:,2) + b(3).*(1+x(:,1)))*1e6+to_fit_val_offset;%full model for how the tune out behaves

surf_samp_to=surf_mdl(fit_mdl_full.Coefficients.Estimate,[reshape(surf_polz_q,[],1),reshape(surf_polz_v,[],1)]);
surf_samp_to=reshape(surf_samp_to,size(surf_polz_q,1),size(surf_polz_q,2));
h=surface(surf_polz_v,surf_polz_q,(surf_samp_to-to_scalar_minus_half_tensor.val)*1e-6);
colormap(viridis)
set(h, 'FaceAlpha', 0.7)
shading interp
hold on
scatter3(polz_v,polz_q,(to_val_for_polz-to_scalar_minus_half_tensor.val)*1e-6,'MarkerFaceColor',[0 .75 .75])
scatter3(0,-1,0*1e-6,'MarkerFaceColor','r')
grid on
title('Full stokes space representation')
xlabel('Fourth Stokes Parameter, V')
ylabel('Second Stokes Parameter, Q')
zlabel( sprintf('Tune out value-%.1f±%.1f (MHz)',to_scalar_minus_half_tensor.val*1e-6,to_scalar_minus_half_tensor.unc*1e-6))
set(gcf,'color','w')
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1600, 900])
zlim([-8,10]*1e3)
view(-25,45)

stfig('2d TO residuals');
to_mdl_val=surf_mdl(fit_mdl_full.Coefficients.Estimate,[polz_q,polz_v]);
to_fit_resid=to_val_for_polz-to_mdl_val;
scatter3(polz_v,polz_q,to_fit_resid*1e-6,'MarkerFaceColor',[0 .75 .75])
hold on
h=surface(surf_polz_v,surf_polz_q,surf_polz_q.*0);
colormap(viridis)
set(h, 'FaceAlpha', 0.7)
shading interp
hold off
zlim([-10,10]*1e3)
view(-25,45)
xlabel('Fourth Stokes Parameter, V')
ylabel('Second Stokes Parameter, Q')
zlabel('Measurments-model (MHz)')

%% Plots of tune-out dependance on the individual stokes parameters
modelfun = @(b,x) b(1) + b(2).*x(:,1);

vec_corr_to = to_val_for_polz-polz_v.*fit_vals(2);
beta0 = [fit_vals(1),fit_vals(3)];
% fit_mdl_t = fitnlm(polz_q,vec_corr_to./1e6,modelfun,beta0,...
%     'Options',opts,'Weights',wlin,'CoefficientNames' ,{'tune_out','tensor'});
fit_mdl_t=full_mdl

tens_corr_to = to_val_for_polz-polz_q.*fit_vals(3);
beta0 = [fit_vals(1),fit_vals(2)];
fit_mdl_v = fitnlm(polz_v,tens_corr_to./1e6,modelfun,beta0,...
    'Options',opts,'Weights',wlin,'CoefficientNames' ,{'tune_out','vector'});

disp_config.fig_name = 'TO fit stokes 4';
disp_config.plot_title='';
disp_config.x_label = 'Fourth Stokes Parameter, V';
disp_config.plot_offset.val=fit_vals(1)./1e6;
disp_config.plot_offset.unc=fit_uncs(1)./1e6;
%plot_binned_nice(disp_config,x_dat,y_dat,weights,fit_mdl)
plot_binned_nice(disp_config,polz_v,tens_corr_to./1e6,wlin,fit_mdl_v)

disp_config.fig_name = 'TO fit stokes 2';
disp_config.x_label = 'Second Stokes Parameter, Q';
disp_config.plot_offset.val=fit_vals(1)./1e6;
disp_config.plot_offset.unc=fit_uncs(1)./1e6;
plot_binned_nice(disp_config,polz_q,vec_corr_to./1e6,wlin,fit_mdl_t)

%% Residuals
to_res_run = to_val_for_polz-predict(fit_mdl,[polz_q,polz_v]);
to_res = to_val_for_polz-predict(fit_mdl,[polz_q,polz_v]);
stfig('residuals')
sfigure(33099);
plot(to_res_run./1e6)
xlabel('index')
ylabel('residual (Mhz)')
sfigure(4445005);
histfit(to_res./1e6,50)
xlabel('residual (MHz)')
ylabel('count')


%%
%to_val_run = [data.hwp.main.lin_fit{1};data.qwp.main.lin_fit{1};data.mix.main.lin_fit{1}];%the to val for a particular run

%V = zeros(numel(to_val),1);%fourth stokes parameter
%d_p = zeros(numel(to_val),1);%power diffrence between min and max transmition
%theta = zeros(numel(to_val),1);%min angle

%[V_run,d_p_run,theta_run] = pol_data_import(pol_opts);

% shot_idx = 1;
% scan_num_vec = [data.hwp.main.scan_num;data.qwp.main.scan_num;data.mix.main.scan_num];
% 
% for loop_idx=1:size(V_run,1)
%     scan_num = scan_num_vec(loop_idx);
%     V(shot_idx:(shot_idx+scan_num-1),1) = ones(scan_num,1).*V_run(loop_idx);
%     d_p(shot_idx:(shot_idx+scan_num-1),1) = ones(scan_num,1).*d_p_run(loop_idx);
%     theta(shot_idx:(shot_idx+scan_num-1),1) = ones(scan_num,1).*theta_run(loop_idx);
%     shot_idx = shot_idx+scan_num;
% end


function to_scalar_minus_half_tensor=two_stage_two_dim_to_fit(full_data)

full_data=cell2mat(full_data);

polz_cont=full_data(:,1);
polz_v=full_data(:,2);
polz_theta=full_data(:,3);
to_val_for_polz=full_data(:,4);
to_unc=full_data(:,5);
wlin=1./(to_unc.^2);




Q_fun = @(d_p,theta,phi) d_p.*cos(2.*(theta+phi.*pi));%function to calculate second stokes parameter, in a rotated frame
full_mdl = @(b,x) b(1) + b(2).*x(:,2) + b(3).*(1+Q_fun(x(:,1),x(:,3),b(4)));%full model for how the tune out behaves
vars = [polz_cont,polz_v,polz_theta];
to_val_fit = to_val_for_polz;
wlin_fit = wlin./sum(wlin);
to_fit_val_offset=wmean(to_val_fit,wlin_fit);

 gf_opt=[];
gf_opt.domain=[[-1,1]*1e5;...   %offset
               [-1,1]*1e4;...  %vector
               [-1,1]*1e4;...  %tensor
               [-1,1]*4*pi;...  %phase
               ];        
gf_opt.start=[0, 1e3, 1*1e3, (rand(1)-0.5)*2*pi];
gf_opt.rmse_thresh=500;
gf_opt.plot=false;
gf_opt.level=2;
gf_out=global_fit(vars,(to_val_fit-to_fit_val_offset)*1e-6,full_mdl,gf_opt);


%
opts = statset('MaxIter',1e4);
beta0 = gf_out.params;%0.2 or 0.8
fit_mdl_full = fitnlm(vars,(to_val_fit-to_fit_val_offset)*1e-6,full_mdl,beta0,...
    'Options',opts,'Weights',wlin_fit,'CoefficientNames' ,{'tune_out','vector','tensor','phase'});
fit_vals_full = fit_mdl_full.Coefficients.Estimate;
%
%second fit with phase fixed
polz_q = Q_fun(polz_cont,polz_theta,fit_vals_full(4));%second stokes parameter, relative to the vector that is otrhogonal to the B feild and beam axis
% TO CHECK
%Q_run = Q_fun(polz_cont,polz_theta,fit_vals(4));%second stokes parameter for each run
full_mdl = @(b,x) b(1) + b(2).*x(:,2) + b(3).*(1+x(:,1));%full model for how the tune out behaves
vars = [polz_q,polz_v];
beta0 = fit_vals_full(1:3);
opts = statset('MaxIter',1e4);
fit_mdl = fitnlm(vars,(to_val_fit-to_fit_val_offset)*1e-6,full_mdl,beta0,...
    'Options',opts,'Weights',wlin_fit,'CoefficientNames' ,{'tune_out','vector','tensor'});
fit_coefs_fixed_phase.val = fit_mdl.Coefficients.Estimate;
fit_coefs_fixed_phase.unc = fit_mdl.Coefficients.SE;

to_scalar_minus_half_tensor=[];
to_scalar_minus_half_tensor=(fit_coefs_fixed_phase.val(1)-0*fit_coefs_fixed_phase.val(3))*1e6+to_fit_val_offset;

end
