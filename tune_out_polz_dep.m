clear all


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
%% import data
% %% HWP
% % loop_config.dir = {
% %     '..\scratch_data\20190227_qwp_270',
% %     '..\scratch_data\20190227_qwp_286',
% %     '..\scratch_data\20190227_qwp_310',
% %     };
% %root_data_dir='..\scratch_data';
% root_data_dir='E:\scratch\good_data';
% files = dir(root_data_dir);
% files=files(3:end);
% % Get a logical vector that tells which is a directory.
% dir_mask = [files.isdir];
% folders=files(dir_mask);
% folders=arrayfun(@(x) fullfile(root_data_dir,x.name),folders,'UniformOutput' ,false);
% folders
% loop_config.dir=folders;
% %%
% 
% data = load_pocessed_to_data(loop_config);

%%
%save('./data/20190611_imported_data_for_tune_out_polz_dep.mat')
load('./data/20190611_imported_data_for_tune_out_polz_dep.mat')

%% Theory
% TO val from https://journals.aps.org/pra/abstract/10.1103/PhysRevA.93.052516, 413.0859(4)
% half tensor shift from \theta_p=54.74 \lambda_TO-> 413.083876e-9
% half tensor shift from \theta_p=90 \lambda_TO-> 413.082896e-9
to_theory=[];
to_theory.wl.val = 413.0859e-9 %+ (413.082896e-9 -413.083876e-9);
to_theory.wl.unc = 0.0004e-9;
[to_theory.freq.val,to_theory.freq.unc] = f2wl(to_theory.wl.val,to_theory.wl.unc); %best theory guess in Hz
fprintf('theory val freq %.0f±%.0f \n',to_theory.freq.val*1e-6,to_theory.freq.unc*1e-6)

to_old.wl.val=413.0938e-9;
to_old.wl.unc=sqrt(0.0009e-9^2+0.0020e-9^2);
[to_old.freq.val,to_old.freq.unc] = f2wl(to_old.wl.val,to_old.wl.unc); 

%% Generate data vectors
%to_val = [data.hwp.drift.to_val{1};data.qwp.drift.to_val{1};data.mix.drift.to_val{1}];%all our single scan tune-out values
%to_unc = [data.hwp.drift.to_val{2};data.qwp.drift.to_val{2};data.mix.drift.to_val{2}];%all corresponding uncertainties
to_val_for_polz=data.drift.to.val;
to_unc=data.drift.to.unc;
wlin=1./(to_unc.^2);

%% polarisation model/data options
pol_opts.location = 'post';%post, pre_cen, pre_left, pre_right
pol_opts.predict_method = 'full_fit_pref_fit';%'full_fit_pref_fit','full_fit_pref_data','full_fit_only','only_data'; %obs (obsovation) fit (pertial fit) full_fit (fit with all parameters free)
                                        %'interp_only','interp_pref_data','interp_pref_interp'
                                        %'gauss_only','gauss_pref_data','gauss_pref_interp'
pol_opts.smoothing=3; %deg
pol_opts.wrap_hwp=0;

pol_opts.hwp=data.drift.wp.hwp;
pol_opts.qwp=data.drift.wp.qwp;
pol_model=pol_data_query(pol_opts);

polz_theta=pol_model.theta.val;
polz_v=pol_model.v.val;
polz_cont=pol_model.cont.val;

to_val_fit = to_val_for_polz;
wlin_fit = wlin./sum(wlin);
to_fit_val_offset=wmean(to_val_fit,wlin_fit);


% set up for the bootstrap wrapper
full_data_set=[polz_cont,polz_v,polz_theta,to_val_for_polz,to_unc];
full_data_set=num2cell(full_data_set,2);
% check that function gives the same answer as the above procedure
[wrap_fun_ans,detailed_out]=two_stage_two_dim_to_fit(full_data_set,to_fit_val_offset,1);
wrap_fun_ans=wrap_fun_ans*1e6+to_fit_val_offset;

to_scalar_minus_half_tensor=detailed_out.tsmht;
fit_mdl_full=detailed_out.fit;
fit_vals_full=fit_mdl_full.Coefficients.Estimate;

fprintf('to val model predict %.1f\n',(to_scalar_minus_half_tensor.val_predict*1e-6-725730000))
fprintf('to val model vals    %.1f\n',(to_scalar_minus_half_tensor.val*1e-6-725730000))
fprintf('to val wav    %.6f\n',f2wl(to_scalar_minus_half_tensor.val)*1e9)
fprintf('to unc model predict %.1f\n',(to_scalar_minus_half_tensor.unc_predict*1e-6))
fprintf('to unc model vals    %.1f\n',(to_scalar_minus_half_tensor.unc*1e-6))
fprintf('theta_k %.1f\n',mod(fit_vals_full(5),pi/2))
fprintf('angle offset %.3f\n',mod(fit_vals_full(4),pi))
fprintf('fit rmse %f\n',fit_mdl_full.RMSE)


%check that the tune out value is recovered with the correct query to the model function
fit_mdl_err=full_tune_out_polz_model(fit_vals_full,[-1,0,-fit_vals_full(4)])*1e6+to_fit_val_offset-to_scalar_minus_half_tensor.val;
if fit_mdl_err>eps
    error('cannot replicate prediction from wrapper function')
end



% Do the bootstrap

% %detailed bootstrap
boot=bootstrap_se(@two_stage_two_dim_to_fit,full_data_set,...
    'opp_arguments',{to_fit_val_offset,1},...
    'plots',true,...
    'replace',true,...
    'plot_fig_name','TO boot MHz',...
    'samp_frac_lims',[0.01,1],...%[0.005,0.9]
    'num_samp_frac',5,... %20
    'num_samp_rep',100,... %1e2
    'plot_fig_name','TO fit bootstrap',...
    'save_multi_out',0,...
    'verbose',3);

% boot=bootstrap_se(@two_stage_two_dim_to_fit,full_data_set,...
%     'opp_arguments',{to_fit_val_offset,1},...
%     'plots',true,...
%     'replace',true,...
%     'plot_fig_name','TO boot MHz',...
%     'samp_frac_lims',1,...%[0.005,0.9]
%     'num_samp_frac',1,... %20
%     'num_samp_rep',200,... %1e2
%     'plot_fig_name','TO fit bootstrap',...
%     'save_multi_out',0,...
%     'verbose',3);

to_scalar_minus_half_tensor.unc_boot=boot.results.se_fun_whole*1e6;
to_scalar_minus_half_tensor.unc_unc_boot=boot.results.se_se_fun_whole*1e6;



%% spit out the results

to_wav_val=f2wl(to_scalar_minus_half_tensor.val);
to_wav_unc=2*to_scalar_minus_half_tensor.unc*const.c/(to_scalar_minus_half_tensor.val^2);
fprintf('\n====TO full fit results==========\n')
fprintf('TO freq             %.1f±(%.0f(fit vals),%.0f(fit predict),%.0f±%.0f(boot)) MHz\n',...
    to_scalar_minus_half_tensor.val*1e-6,...
    to_scalar_minus_half_tensor.unc*1e-6,...
    to_scalar_minus_half_tensor.unc_predict*1e-6,...
    to_scalar_minus_half_tensor.unc_boot*1e-6,...
    to_scalar_minus_half_tensor.unc_unc_boot*1e-6)
fprintf('diff from TOV1      %.0f±%.0f MHz \n',(to_scalar_minus_half_tensor.val-to_old.freq.val)*1e-6,sqrt(to_scalar_minus_half_tensor.unc^2+to_old.freq.unc^2)*1e-6)
fprintf('diff from Theory    %.0f±%.0f MHz \n',...
    (to_scalar_minus_half_tensor.val-to_theory.freq.val)*1e-6,...
    sqrt(to_scalar_minus_half_tensor.unc^2+to_theory.freq.unc^2)*1e-6)
fprintf('TO wavelength       %.6f±%f nm \n',to_wav_val*1e9,to_wav_unc*1e9)

%% plot options

font_name='cmr10';
font_size_global=20;
folt_size_label=20;


%%
% Full 2D plot of the fit in 2d stokes space, with each scan shown
% we will plot in the second (Q) and fourth (v) stokes parameters, 
% we rotate the measured stokes parameters into the convinent basis found from the above fit

%we will compute the second stokes parameter for the data
polz_q=Q_fun(polz_cont,polz_theta,fit_vals_full(4));

stfig('2d TO fit');
clf
[surf_polz_q, surf_polz_v] = meshgrid(linspace(-1,1),linspace(-1,1));
surf_mdl_fun = @(b,x) b(1) + (1/2).*x(:,2).*cos(b(5)).*b(2) - (1/2)*D_fun(b(5),x(:,1)).*b(3);%full model for how the tune out behaves
surf_samp_to=surf_mdl_fun(fit_mdl_full.Coefficients.Estimate,[reshape(surf_polz_q,[],1),reshape(surf_polz_v,[],1)]);
surf_samp_to=surf_samp_to*1e6+to_fit_val_offset;
surf_samp_to=reshape(surf_samp_to,size(surf_polz_q,1),size(surf_polz_q,2));
h=surface(surf_polz_v,surf_polz_q,(surf_samp_to-to_scalar_minus_half_tensor.val)*1e-6);
colormap(viridis)
set(h, 'FaceAlpha', 0.7)
shading interp
hold on
scatter3(0,-1,0*1e-6,'MarkerFaceColor','r')

%check that we get the same result using the fit model directly, and use this for CI
% keep in mind what we fit looks like this
%result=b(1) + (1/2).*x(:,2).*cos(b(5)).*b(2) - (1/2)*D_fun(b(5),Q_fun(x(:,1),x(:,3),b(4))).*b(3);
% where b=tune_out_scalar,reduced_vector,reduce_tensor,angle between polz measurment basis and B cross k, theta k
% and x=polz_cont,polz_v,polz_theta
% where Q fun
% Q_fun(d_p,theta,phi)=d_p.*cos(2.*(theta+phi));
% i want to pass x to evaluate the model for some target value of V,Q, a slightly hacky approach
% if we pass in -1*fit_vals_full(4) as 3rd precictor which
% becomes the 2nd input of Q_fun and will give d_p*cos(0),the first argument of Q_fun is the value it returns
predict_predictor_input=cat(2,reshape(surf_polz_q,[],1),reshape(surf_polz_v,[],1),repmat(-fit_vals_full(4),numel(surf_polz_v),1));
[surf_mdl_val,surf_mdl_ci]=predict(fit_mdl_full,predict_predictor_input,'Alpha',1-erf(1/sqrt(2)),'Prediction','Curve');
%'Prediction','Curve','Prediction','observation'
surf_mdl_val=surf_mdl_val*1e6+to_fit_val_offset;
surf_mdl_val=reshape(surf_mdl_val,size(surf_polz_q,1),size(surf_polz_q,2));

if max(reshape(surf_mdl_val-surf_samp_to,[],1))>eps(1)
    error('model values from prediction and manual method do not agree')
end
% hold on
% h=surface(surf_polz_v,surf_polz_q,(surf_mdl_val-to_scalar_minus_half_tensor.val)*1e-6);
% set(h, 'FaceAlpha', 0.3)
% shading interp


surf_mdl_ci=surf_mdl_ci*1e6+to_fit_val_offset;
surf_mdl_ci=reshape(surf_mdl_ci,size(surf_polz_q,1),size(surf_polz_q,2),2);

h=surface(surf_polz_v,surf_polz_q,(surf_mdl_ci(:,:,1)-to_scalar_minus_half_tensor.val)*1e-6);
colormap(viridis)
set(h, 'FaceAlpha', 0.7)
shading interp

h=surface(surf_polz_v,surf_polz_q,(surf_mdl_ci(:,:,2)-to_scalar_minus_half_tensor.val)*1e-6);
colormap(viridis)
set(h, 'FaceAlpha', 0.7)
shading interp


scatter3(polz_v,polz_q,(to_val_for_polz-to_scalar_minus_half_tensor.val)*1e-6,'MarkerFaceColor',[0 .75 .75])
grid on
title('Full stokes space representation')
xlabel('Fourth Stokes Parameter, V')
ylabel('Second Stokes Parameter, Q')
zlabel( sprintf('Tune out value-%.1f±%.1f (MHz)',to_scalar_minus_half_tensor.val*1e-6,to_scalar_minus_half_tensor.unc_boot*1e-6))
set(gcf,'color','w')
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1600, 900])
zlim([-8,10]*1e3)
view(-25,45)

stfig('2d TO residuals');
to_mdl_val=surf_mdl_fun(fit_mdl_full.Coefficients.Estimate,[polz_q,polz_v]);
to_mdl_val=to_mdl_val*1e6+to_fit_val_offset;
to_fit_resid=to_val_for_polz-to_mdl_val;
scatter3(polz_v,polz_q,to_fit_resid*1e-6,'MarkerFaceColor',[0 .75 .75])
hold on
h=surface(surf_polz_v,surf_polz_q,surf_polz_q.*0);
colormap(viridis)
set(h, 'FaceAlpha', 0.7)
shading interp
hold off
zlim([-4,4]*1e3)
view(-25,45)
xlabel('Fourth Stokes Parameter, V')
ylabel('Second Stokes Parameter, Q')
zlabel('Measurments-model (MHz)')

%% Full 2D plot of the fit in 2d stokes space, with each run binned into a single point
% we will plot in the second (Q) and fourth (v) stokes parameters


colors_main=[[53,126,220];[33,188,44];[0,0,0]]./255;
lch=colorspace('RGB->LCH',colors_main(:,:));
lch(:,1)=lch(:,1)+20;
colors_detail=colorspace('LCH->RGB',lch);
%would prefer to use srgb_2_Jab here
color_shaded=colorspace('RGB->LCH',colors_main(3,:));
color_shaded(1)=50;
color_shaded=colorspace('LCH->RGB',color_shaded);


stfig('2d TO fit, binned');
clf
surf_mdl=surface(surf_polz_v,surf_polz_q,(surf_samp_to-to_scalar_minus_half_tensor.val)*1e-6);
set(surf_mdl,'linestyle','none')
%set(surf_mdl,'FaceColor',[255 127 42]./255)
colormap(viridis)
set(surf_mdl, 'FaceAlpha', 0.5)
hold on
% plot the CI surfaces
% surf_ci=surface(surf_polz_v,surf_polz_q,(surf_mdl_ci(:,:,1)-to_scalar_minus_half_tensor.val)*1e-6);
% set(surf_ci,'linestyle','none')
% set(surf_ci,'FaceColor',[0 0 0])
% set(surf_ci, 'FaceAlpha', 0.3)
% surf_ci=surface(surf_polz_v,surf_polz_q,(surf_mdl_ci(:,:,2)-to_scalar_minus_half_tensor.val)*1e-6);
% set(surf_ci,'linestyle','none')
% set(surf_ci,'FaceColor',[0 0 0])
% set(surf_ci, 'FaceAlpha', 0.3)

%bin the measured data
polz_vq=[polz_v,polz_q];
unique_polz_vq = unique(polz_vq, 'rows');
unique_polz_vq=unique_polz_vq(sum(isnan(unique_polz_vq),2)==0,:);

% ok this line is very dense
% use col_row_fun_mat to run the operation on rows of the matrix unique_polz_vq
% i select the elements of to_val_for_polz that corespond to matching the polz_vq matrix rows and then
% compute the unc_wmean which returns a vector of weighted mean,weighted unc & unweighted std
polz_vq_to_val_unc=col_row_fun_mat(@(x) unc_wmean_vec(to_val_for_polz(all(x'==polz_vq,2)),to_unc(all(x'==polz_vq,2))),...
                unique_polz_vq,2);
polz_vq_to_val_unc(:,3:5)=polz_vq_to_val_unc;         
polz_vq_to_val_unc(:,1:2)=unique_polz_vq; 
%scale
polz_vq_to_val_unc_scaled=polz_vq_to_val_unc;
polz_vq_to_val_unc_scaled(:,3)=(polz_vq_to_val_unc(:,3)-to_scalar_minus_half_tensor.val)*1e-6;
polz_vq_to_val_unc_scaled(:,4)=polz_vq_to_val_unc(:,4)*1e-6; %scale the unc
polz_vq_to_val_unc_scaled(:,5)=polz_vq_to_val_unc(:,5)*1e-6; %scale the std

plot3d_errorbars(polz_vq_to_val_unc_scaled(:,1), polz_vq_to_val_unc_scaled(:,2), polz_vq_to_val_unc_scaled(:,3), [], [], polz_vq_to_val_unc_scaled(:,4),colors_main(1,:));

scatter3(polz_vq_to_val_unc_scaled(:,1), polz_vq_to_val_unc_scaled(:,2), polz_vq_to_val_unc_scaled(:,3),50,'o',...
   'MarkerEdgeColor',colors_main(1,:),'MarkerFaceColor',colors_detail(1,:),'LineWidth',2)
scatter3(0,-1,0,150,'s','MarkerFaceColor','k','MarkerEdgeColor','k')
grid on
%title('Full stokes space representation')
xlabel('\bfv\rm')
ylabel('\bfq\rm')
%zlabel( sprintf('Tune out value-%.1f±%.1f (MHz)',to_scalar_minus_half_tensor.val*1e-6,to_scalar_minus_half_tensor.unc_boot*1e-6))
zlabel('\omega_{TO}(\bfq\rm,\bfv\rm) - \omega^{SMHT}_{TO} (MHz)')
set(gcf,'color','w')
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1600, 900])
zlim([-8,10]*1e3)
view(-25,45)
hold off
axis vis3d

stfig('2d TO residuals, binned');
polz_vq_to_resid_unc=polz_vq_to_val_unc;
to_mdl_val=surf_mdl_fun(fit_mdl_full.Coefficients.Estimate,polz_vq_to_resid_unc(:,[2,1]));
to_mdl_val=to_mdl_val*1e6+to_fit_val_offset;

polz_vq_to_resid_unc(:,3)=polz_vq_to_resid_unc(:,3)-to_mdl_val;

polz_vq_to_resid_unc_scaled=polz_vq_to_resid_unc;
polz_vq_to_resid_unc_scaled(:,[3,4])=polz_vq_to_resid_unc_scaled(:,[3,4])*1e-6;

scatter3(polz_vq_to_resid_unc_scaled(:,1), polz_vq_to_resid_unc_scaled(:,2), polz_vq_to_resid_unc_scaled(:,3),'MarkerFaceColor',[0 .75 .75])
hold on
plot3d_errorbars(polz_vq_to_resid_unc_scaled(:,1), polz_vq_to_resid_unc_scaled(:,2), polz_vq_to_resid_unc_scaled(:,3), [], [], polz_vq_to_resid_unc_scaled(:,4),'k');
surf_mdl=surface(surf_polz_v,surf_polz_q,surf_polz_q.*0);
set(surf_mdl,'linestyle','none')
set(surf_mdl,'FaceColor',[255 127 42]./255)
set(surf_mdl, 'FaceAlpha', 0.5)

surf_ci=surface(surf_polz_v,surf_polz_q,(surf_mdl_ci(:,:,1)-surf_samp_to)*1e-6);
set(surf_ci,'linestyle','none')
set(surf_ci,'FaceColor',[0 0 0])
set(surf_ci, 'FaceAlpha', 0.3)
surf_ci=surface(surf_polz_v,surf_polz_q,(surf_mdl_ci(:,:,2)-surf_samp_to)*1e-6);
set(surf_ci,'linestyle','none')
set(surf_ci,'FaceColor',[0 0 0])
set(surf_ci, 'FaceAlpha', 0.3)


hold off
zlim([-1,1]*1e3)
view(-25,45)
xlabel('Fourth Stokes Parameter, V')
ylabel('Second Stokes Parameter, Q')
zlabel('Measurments-model (MHz)')



%% Plots of tune-out dependance on the individual stokes parameters

% plot with Q, stokes 2nd
samp_q=col_vec(linspace(-1.1,1.1,1e4));
predict_predictor_input=cat(2,samp_q,samp_q*0,repmat(-fit_vals_full(4),numel(samp_q),1));
[fit_mdl_val,fit_mdl_ci]=predict(fit_mdl_full,predict_predictor_input,'Alpha',1-erf(1/sqrt(2)),'Prediction','Curve');
%'Prediction','Curve','Prediction','observation'
fit_mdl_val=fit_mdl_val*1e6+to_fit_val_offset;

fit_mdl_ci=fit_mdl_ci*1e6+to_fit_val_offset;
stfig('Q dep');
clf
patch([samp_q', fliplr(samp_q')],...
    ([fit_mdl_ci(:,1)', fliplr(fit_mdl_ci(:,2)')]-to_scalar_minus_half_tensor.val).*1e-6,...
    [1,1,1].*0.7,'EdgeColor','none')
hold on
plot(samp_q,(fit_mdl_val-to_scalar_minus_half_tensor.val)*1e-6,'k')
xlabel('Q value')
ylabel('Tune out -TOSMHT (MHz)')
title('V=0 extrapolation')
% now I want to plot every scan with its error bar that has been corrected onto V=0

shift_vq=v_correcting_shift(polz_v,polz_q,fit_mdl_full);

shifted_scaled_to_vals=(to_val_for_polz-to_scalar_minus_half_tensor.val+shift_vq(:,1)*1e6)*1e-6;
errorbar(polz_q,shifted_scaled_to_vals,to_unc*1e-6...
     ,'o','CapSize',0,'MarkerSize',5,'Color','r',...
     'LineWidth',1.5);
%ylim([min(shifted_scaled_to_vals),max(shifted_scaled_to_vals)])
ylim([-0.5,0.5]*1e4)
xlim([-1.1,1.1])

% do the same for V, stokes 4th
samp_v=col_vec(linspace(-1.1,1.1,1e4));
predict_predictor_input=cat(2,samp_v*0-1,samp_v,repmat(-fit_vals_full(4),numel(samp_v),1));
[fit_mdl_val,fit_mdl_ci]=predict(fit_mdl_full,predict_predictor_input,'Alpha',1-erf(1/sqrt(2)),'Prediction','Curve');
%'Prediction','Curve','Prediction','observation'
fit_mdl_val=fit_mdl_val*1e6+to_fit_val_offset;
fit_mdl_ci=fit_mdl_ci*1e6+to_fit_val_offset;
stfig('V dep');
clf
patch([samp_v', fliplr(samp_v')],...
    ([fit_mdl_ci(:,1)', fliplr(fit_mdl_ci(:,2)')]-to_scalar_minus_half_tensor.val).*1e-6,...
    [1,1,1].*0.7,'EdgeColor','none')
hold on
plot(samp_v,(fit_mdl_val-to_scalar_minus_half_tensor.val)*1e-6,'k')
xlabel('V value')
ylabel('Tune out -TOSMHT(MHz)')
title('Q=-1 extrapolation')
shifted_scaled_to_vals=(to_val_for_polz-to_scalar_minus_half_tensor.val+shift_vq(:,2)*1e6)*1e-6;

errorbar(polz_v,shifted_scaled_to_vals,to_unc*1e-6...
     ,'o','CapSize',0,'MarkerSize',5,'Color','r',...
     'LineWidth',1.5);
%ylim([min(shifted_scaled_to_vals),max(shifted_scaled_to_vals)]) 
ylim([-1,0.8].*1e4) 
xlim([-1.1,1.1])
 

%% Plot the V,Q plots with binned data

colors_main=[[53,126,220];[33,188,44];[0,0,0]]./255;
lch=colorspace('RGB->LCH',colors_main(:,:));
lch(:,1)=lch(:,1)+20;
colors_detail=colorspace('LCH->RGB',lch);
%would prefer to use srgb_2_Jab here
color_shaded=colorspace('RGB->LCH',colors_main(3,:));
color_shaded(1)=50;
color_shaded=colorspace('LCH->RGB',color_shaded);


% plot with Q, stokes 2nd
samp_q=col_vec(linspace(-1.1,1.1,1e4));
predict_predictor_input=cat(2,samp_q,samp_q*0,repmat(-fit_vals_full(4),numel(samp_q),1));
[fit_mdl_val,fit_mdl_ci]=predict(fit_mdl_full,predict_predictor_input,'Alpha',1-erf(1/sqrt(2)),'Prediction','Curve');
%'Prediction','Curve','Prediction','observation'
fit_mdl_val=fit_mdl_val*1e6+to_fit_val_offset;

fit_mdl_ci=fit_mdl_ci*1e6+to_fit_val_offset;
stfig('Q dep,binned');
clf
patch([samp_q', fliplr(samp_q')],...
    ([fit_mdl_ci(:,1)', fliplr(fit_mdl_ci(:,2)')]-to_scalar_minus_half_tensor.val).*1e-6,...
    [1,1,1].*0.7,'EdgeColor','none')
hold on
plot(samp_q,(fit_mdl_val-to_scalar_minus_half_tensor.val)*1e-6,'k')
xlabel('2^{nd} stokes parameter, \bfq\rm')
ylabel('\omega_{TO}(\bfq\rm,\bfv\rm=0) -\omega^{SMHT}_{TO} (MHz)')
%title('v=0 extrapolation')
% now I want to plot every scan with its error bar that has been corrected onto V=0

shift_vq=v_correcting_shift(polz_vq_to_val_unc(:,1),polz_vq_to_val_unc(:,2),fit_mdl_full);

shifted_scaled_to_vals=(polz_vq_to_val_unc(:,3)-to_scalar_minus_half_tensor.val+shift_vq(:,1)*1e6)*1e-6;


errorbar(polz_vq_to_val_unc(:,2),shifted_scaled_to_vals,polz_vq_to_val_unc(:,5)*1e-6...
     ,'o','CapSize',0,'Marker','none','Color',colors_detail(1,:),...
     'LineWidth',1.5);
errorbar(polz_vq_to_val_unc(:,2),shifted_scaled_to_vals,polz_vq_to_val_unc(:,4)*1e-6,...
    'o','CapSize',0,'MarkerSize',5,'Color',colors_main(1,:),...
    'MarkerFaceColor',colors_detail(1,:),'LineWidth',2.8);
%ylim([min(shifted_scaled_to_vals),max(shifted_scaled_to_vals)])
ylim([-2.5,1.2]*1e3)
xlim([-1.1,1.1])
hold off
set(gca,'FontSize',font_size_global,'FontName',font_name)
set(gcf,'Units','Pixels')
set(gcf,'Position',[1068         355         676         453])




stfig('Q dep,binned,residuals');
subplot(2,1,1)
samp_q=polz_vq_to_val_unc(:,2);
predict_predictor_input=cat(2,samp_q,samp_q*0,repmat(-fit_vals_full(4),numel(samp_q),1));
fit_mdl_val=predict(fit_mdl_full,predict_predictor_input,'Alpha',1-erf(1/sqrt(2)),'Prediction','Curve');
%'Prediction','Curve','Prediction','observation'
fit_mdl_val=fit_mdl_val*1e6+to_fit_val_offset;
shifted_scaled_to_vals=(polz_vq_to_val_unc(:,3)+shift_vq(:,1)*1e6-fit_mdl_val)*1e-6;
errorbar(samp_q,shifted_scaled_to_vals,polz_vq_to_val_unc(:,5)*1e-6...
     ,'o','CapSize',0,'Marker','none','Color','g',...
     'LineWidth',1.5);
hold on
errorbar(samp_q,shifted_scaled_to_vals,polz_vq_to_val_unc(:,4)*1e-6...
     ,'o','CapSize',0,'MarkerSize',5,'Color','r',...
     'LineWidth',1.5);
hold off
subplot(2,1,2)
residuals_num_ste=shifted_scaled_to_vals./(polz_vq_to_val_unc(:,4)*1e-6);
plot(samp_q,residuals_num_ste,'ok')
fprintf('sd of (residuals/ste) in Q dep %.1f \n',std(residuals_num_ste))
xlabel('2^{nd} stokes parameter, Q')
ylabel('Error From Model (MHz)')
set(gca,'FontSize',font_size_global,'FontName',font_name)





% do the same for V, stokes 4th
samp_v=col_vec(linspace(-1.1,1.1,1e4));
predict_predictor_input=cat(2,samp_v*0-1,samp_v,repmat(-fit_vals_full(4),numel(samp_v),1));
[fit_mdl_val,fit_mdl_ci]=predict(fit_mdl_full,predict_predictor_input,'Alpha',1-erf(1/sqrt(2)),'Prediction','Curve');
%'Prediction','Curve','Prediction','observation'
fit_mdl_val=fit_mdl_val*1e6+to_fit_val_offset;
fit_mdl_ci=fit_mdl_ci*1e6+to_fit_val_offset;
stfig('V dep,binned');
clf
patch([samp_v', fliplr(samp_v')],...
    ([fit_mdl_ci(:,1)', fliplr(fit_mdl_ci(:,2)')]-to_scalar_minus_half_tensor.val).*1e-6,...
    [1,1,1].*0.7,'EdgeColor','none')
hold on
plot(samp_v,(fit_mdl_val-to_scalar_minus_half_tensor.val)*1e-6,'k')
xlabel('4^{th} stokes parameter,\bfv\rm')
ylabel('\omega_{TO}(\bfq\rm,=-1,\bfv\rm) -\omega^{SMHT}_{TO} (MHz)')
%title('Q=-1 extrapolation')

shifted_scaled_to_vals=(polz_vq_to_val_unc(:,3)-to_scalar_minus_half_tensor.val+shift_vq(:,2)*1e6)*1e-6;
errorbar(polz_vq_to_val_unc(:,1),shifted_scaled_to_vals,polz_vq_to_val_unc(:,5)*1e-6...
     ,'o','CapSize',0,'Marker','none','Color',colors_detail(1,:),...
     'LineWidth',1.5);
errorbar(polz_vq_to_val_unc(:,1),shifted_scaled_to_vals,polz_vq_to_val_unc(:,4)*1e-6,...
    'o','CapSize',0,'MarkerSize',5,'Color',colors_main(1,:),...
    'MarkerFaceColor',colors_detail(1,:),'LineWidth',2.8);

%ylim([min(shifted_scaled_to_vals),max(shifted_scaled_to_vals)]) 
ylim([-8,8].*1e3) 
xlim([-1.1,1.1])
set(gca,'FontSize',font_size_global,'FontName',font_name)
set(gcf,'Units','Pixels')
set(gcf,'Position',[ 1068         355         676         453])




stfig('V dep,binned,residuals');
subplot(2,1,1)
samp_v=polz_vq_to_val_unc(:,1);
predict_predictor_input=cat(2,samp_v*0-1,samp_v,repmat(-fit_vals_full(4),numel(samp_v),1));
fit_mdl_val=predict(fit_mdl_full,predict_predictor_input,'Alpha',1-erf(1/sqrt(2)),'Prediction','Curve');
%'Prediction','Curve','Prediction','observation'
fit_mdl_val=fit_mdl_val*1e6+to_fit_val_offset;
shifted_scaled_to_vals=(polz_vq_to_val_unc(:,3)+shift_vq(:,2)*1e6-fit_mdl_val)*1e-6;
errorbar(samp_v,shifted_scaled_to_vals,polz_vq_to_val_unc(:,5)*1e-6...
     ,'o','CapSize',0,'Marker','none','Color','g',...
     'LineWidth',1.5);
hold on
errorbar(samp_v,shifted_scaled_to_vals,polz_vq_to_val_unc(:,4)*1e-6...
     ,'o','CapSize',0,'MarkerSize',5,'Color','r',...
     'LineWidth',1.5);
hold off
subplot(2,1,2)
residuals_num_ste=shifted_scaled_to_vals./(polz_vq_to_val_unc(:,4)*1e-6);
plot(polz_vq_to_val_unc(:,2),residuals_num_ste,'ok')
fprintf('sd of (residuals/ste) in V dep %.1f \n',std(residuals_num_ste))



%%

function result=Q_fun(d_p,theta,phi)
result=d_p.*cos(2.*(theta+phi));
%fprintf('Q val %f\n',result);
end

function result=D_fun(theta_k,Q)
result=(3*(sin(theta_k)^2)*((1/2)+(Q/2)))-1;
%fprintf('D val %f\n',result);
end


function result=full_tune_out_polz_model(b,x)
%tune_out_scalar
%reduced_vector
%reduce_tensor
%angle between polz measurment basis and B cross k
% theta k

result=b(1) + (1/2).*x(:,2).*cos(b(5)).*b(2) - (1/2)*D_fun(b(5),Q_fun(x(:,1),x(:,3),b(4))).*b(3);
end


function shift_vq=v_correcting_shift(v_in,q_in,mdl)

q_in=col_vec(q_in);
v_in=col_vec(v_in);
if numel(v_in)~=numel(q_in)
    error('v,q must be same size')
end

predict_predictor_input=cat(2,q_in,v_in,repmat(-mdl.Coefficients.Estimate(4),numel(v_in),1));
fit_mdl_val=predict(mdl,predict_predictor_input);

predict_predictor_input=cat(2,q_in,v_in*0,repmat(-mdl.Coefficients.Estimate(4),numel(v_in),1));
fit_mdl_v_zero=predict(mdl,predict_predictor_input);

predict_predictor_input=cat(2,q_in*0-1,v_in,repmat(-mdl.Coefficients.Estimate(4),numel(v_in),1));
fit_mdl_q_zero=predict(mdl,predict_predictor_input);

shift_vq=[fit_mdl_v_zero-fit_mdl_val,fit_mdl_q_zero-fit_mdl_val];

end

function [to_val_out,detailed_out]=two_stage_two_dim_to_fit(full_data,fit_offset,use_fit_weights)
%return the value in MHZ away from fit_offset as to_val_out
% return more details in detailed_out

full_data=cell2mat(full_data);

polz_cont=full_data(:,1);
polz_v=full_data(:,2);
polz_theta=full_data(:,3);
to_val_for_polz=full_data(:,4);
to_unc=full_data(:,5);


predictor = [polz_cont,polz_v,polz_theta];
to_val_fit = to_val_for_polz;

w_fit=1./(to_unc.^2);
if ~use_fit_weights
    w_fit=w_fit*0+1;
end
w_fit = w_fit./sum(w_fit);

% so here we will use a single peice of information from the atomic theory, namely the gradient of the reduced tensor
% term, from Li-Yan Tang's document "Response to Bryce's questions" \lambda_to (theta_p=0)>\lambda_to (theta_p=90)
% & as the polarizability decreased with increasing wavelength d \alpha^Scalar(\lambda)/ d \lambda <0 
% ? a increase in angle theta_p decreases the polarizability
% & as the prefactor in front of the tensor term decreases with increasing theta_p 
% ? \alpha^Tensor(\omega?\omega_TO)>0
% as \alpha^Scalar(\omega)/ d \omega >0
%  ? the reduced tensor term is positive
tune_out_vals_scaled=(to_val_fit-fit_offset)*1e-6;



gf_opt=[];
gf_opt.domain=[[-1,1]*1e5;...   %tune_out_scalar
               [-1,1]*1e6;...  %reduced_vector*cos_theta_k
               [0,1]*1e6;...  %reduce_tensor   
               [-1,1]*pi/2;...  %angle between polz measurment basis and B cross k
               pi+[-1,1]*pi/4;... % theta k
               ];        
gf_opt.start=[0, 1e3, 1e4, (rand(1)-0.5)*pi/4,pi+(rand(1)-0.5)*pi/4];
gf_opt.rmse_thresh=2e3;
gf_opt.plot=false;
gf_opt.level=2;
gf_out=global_fit(predictor,tune_out_vals_scaled,@full_tune_out_polz_model,gf_opt);

%
opts = statset('MaxIter',1e4);
beta0 = gf_out.params;%0.2 or 0.8
%beta0=[0, 1e3, 1*1e3, (rand(1)-0.5)*2*pi]
fit_mdl_full = fitnlm(predictor,tune_out_vals_scaled,@full_tune_out_polz_model,beta0,...
   'Weights',w_fit,...
    'CoefficientNames' ,{'tune_out_scalar','reduced_vector','reduce_tensor','phase','thetak'},...
    'Options',opts);

fit_vals_full = fit_mdl_full.Coefficients.Estimate;
fit_uncs_full = fit_mdl_full.Coefficients.SE;
%


if fit_vals_full(3)<0
    error('fit reduced tensor term is less than zero')
end


%
[tsmht_details.val_predict,tsmht_details.unc_predict]=...
    predict(fit_mdl_full,[-1,0,-fit_vals_full(4)],'Alpha',1-erf(1/sqrt(2)),'Prediction','Curve');
tsmht_details.val_predict=tsmht_details.val_predict*1e6+fit_offset;
tsmht_details.unc_predict=1/2*range(tsmht_details.unc_predict)*1e6;

tsmht_details.val_no_offset_mhz=fit_vals_full(1)+(1/2)*fit_vals_full(3);
tsmht_details.val=tsmht_details.val_no_offset_mhz*1e6+fit_offset;
tsmht_details.unc=sqrt(fit_uncs_full(1)^2 + (fit_uncs_full(3)/2)^2);
tsmht_details.unc=tsmht_details.unc*1e6;

detailed_out.tsmht=tsmht_details;
detailed_out.fit=fit_mdl_full;

if fit_vals_full(3)<0
    warning('fit has returned a reduced tensor term is less than zero')
    to_val_out=nan;
else
    to_val_out=tsmht_details.val_no_offset_mhz; %return this value in mhz for the bootstrap
end

end








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

% %%
% %second fit with phase fixed
% polz_q = Q_fun(polz_cont,polz_theta,fit_vals_full(4));%second stokes parameter, relative to the vector that is otrhogonal to the B feild and beam axis
% % TO CHECK
% %Q_run = Q_fun(polz_cont,polz_theta,fit_vals(4));%second stokes parameter for each run
% fixed_phase_model = @(b,x) b(1) + (1/2).*x(:,2).*cos(fit_vals_full(4)).*b(2) - (1/2)*D_fun(b(4),Q_fun(x(:,1),x(:,3),fit_vals_full(4))).*b(3);%full model for how the tune out behaves
% predictor = [polz_q,polz_v];
% beta0 = fit_vals_full(1:3);
% opts = statset('MaxIter',1e4);
% fit_mdl = fitnlm(predictor,(to_val_fit-to_fit_val_offset)*1e-6,fixed_phase_model,beta0,...
%     'Options',opts,'Weights',wlin_fit,'CoefficientNames' ,{'tune_out_scalar','reduced_vector','reduce_tensor'})
% fit_coefs_fixed_phase.val = fit_mdl.Coefficients.Estimate;
% fit_coefs_fixed_phase.unc = fit_mdl.Coefficients.SE;
% 
% to_scalar_minus_half_tensor=[];
% to_scalar_minus_half_tensor.val=(fit_coefs_fixed_phase.val(1)-0*fit_coefs_fixed_phase.val(3))*1e6+to_fit_val_offset;
% to_scalar_minus_half_tensor.unc=fit_coefs_fixed_phase.unc(1)*1e6;
% %to_scalar_minus_half_tensor.val=(fit_coefs_fixed_phase.val(1)+1*fit_coefs_fixed_phase.val(3))*1e6+to_fit_val_offset;
% %to_scalar_minus_half_tensor.unc=sqrt(fit_coefs_fixed_phase.unc(1)^2+(0.5*fit_coefs_fixed_phase.unc(3))^2)*1e6;
% fprintf('%.1f\n',(to_scalar_minus_half_tensor.val*1e-6-725735000))
% 
% 
% % 2   2000
% % 0.5 1000
% % 2/3 700
% % Write out the results


% %% OLD WAY
% modelfun = @(b,x) b(1) + b(2).*x(:,1);
% 
% vec_corr_to = to_val_for_polz-polz_v.*fit_vals(2);
% beta0 = [fit_vals(1),fit_vals(3)];
% % fit_mdl_t = fitnlm(polz_q,vec_corr_to./1e6,modelfun,beta0,...
% %     'Options',opts,'Weights',wlin,'CoefficientNames' ,{'tune_out','tensor'});
% fit_mdl_t=full_mdl
% 
% tens_corr_to = to_val_for_polz-polz_q.*fit_vals(3);
% beta0 = [fit_vals(1),fit_vals(2)];
% fit_mdl_v = fitnlm(polz_v,tens_corr_to./1e6,modelfun,beta0,...
%     'Options',opts,'Weights',wlin,'CoefficientNames' ,{'tune_out','vector'});
% 
% disp_config.fig_name = 'TO fit stokes 4';
% disp_config.plot_title='';
% disp_config.x_label = 'Fourth Stokes Parameter, V';
% disp_config.plot_offset.val=fit_vals(1)./1e6;
% disp_config.plot_offset.unc=fit_uncs(1)./1e6;
% %plot_binned_nice(disp_config,x_dat,y_dat,weights,fit_mdl)
% plot_binned_nice(disp_config,polz_v,tens_corr_to./1e6,wlin,fit_mdl_v)
% 
% disp_config.fig_name = 'TO fit stokes 2';
% disp_config.x_label = 'Second Stokes Parameter, Q';
% disp_config.plot_offset.val=fit_vals(1)./1e6;
% disp_config.plot_offset.unc=fit_uncs(1)./1e6;
% plot_binned_nice(disp_config,polz_q,vec_corr_to./1e6,wlin,fit_mdl_t)
% 
% %% Residuals
% to_res_run = to_val_for_polz-predict(fit_mdl,[polz_q,polz_v]);
% to_res = to_val_for_polz-predict(fit_mdl,[polz_q,polz_v]);
% stfig('residuals')
% sfigure(33099);
% plot(to_res_run./1e6)
% xlabel('index')
% ylabel('residual (Mhz)')
% sfigure(4445005);
% histfit(to_res./1e6,50)
% xlabel('residual (MHz)')
% ylabel('count')