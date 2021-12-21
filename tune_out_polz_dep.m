%clear all


%% Exp parameters
aom_freq =0; %aom offset in Hz, as of 20190528 offset is calculated in main_trap_freq

%%
addpath('./lib/Core_BEC_Analysis/lib/') %add the path to set_up_project_path
set_up_project_path('.')


%%
hebec_constants %call the constants function that makes some globals


% import data
%% HWP
% loop_config.dir = {
%     '..\scratch_data\20190227_qwp_270',
%     '..\scratch_data\20190227_qwp_286',
%     '..\scratch_data\20190227_qwp_310',
%     };
%root_data_dir='..\scratch_data';
%root_data_dir='Z:\DataBackup\Bryce_Data_Backup\TO_working_data\to_main_data';
root_data_dir='Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\to_main_data'
files = dir(root_data_dir);
files=files(3:end);
% Get a logical vector that tells which is a directory.
dir_mask = [files.isdir];
folders=files(dir_mask);
folders=arrayfun(@(x) fullfile(root_data_dir,x.name),folders,'UniformOutput' ,false);
folders
loop_config.dir=folders;
%%

%data = load_pocessed_to_data(loop_config);


%%
%save('./data/20211207_imported_data_for_tune_out_polz_dep.mat','data')
load('./data/20211207_imported_data_for_tune_out_polz_dep.mat','data')

%load('./data/20190611_imported_data_for_tune_out_polz_dep.mat')
%load('./data/20211221_imported_data_for_tune_out_polz_dep.mat','data')

%% dist of atom number
% qe=0.09;
% [mean_atomnum,std_atomnum]=pooled_mean_and_std(data.drift.atom_num_probe.val,data.drift.atom_num_probe.std,data.drift.probe_shots);
% mean_atomnum=mean_atomnum/qe;
% std_atomnum=std_atomnum/qe;
% 
% smooth_hist(data.drift.atom_num_probe.val/qe,'sigma',1.5e4)
% hold on
% smooth_hist(data.drift.atom_num_probe.val/qe,'sigma',3e3)
% hold off
% xlabel('atom number')
% ylabel('probability')

%%


%% Theory
% TO val from https://journals.aps.org/pra/abstract/10.1103/PhysRevA.93.052516, 413.0859(4)
% half tensor shift from \theta_p=54.74 \lambda_TO-> 413.083876e-9
% half tensor shift from \theta_p=90 \lambda_TO-> 413.082896e-9
to_theory=[];
to_theory.freq.val=725736430e6;
to_theory.freq.unc=50;
[to_theory.wl.val,to_theory.wl.unc] = f2wl(to_theory.freq.val,to_theory.freq.unc); %best theory guess in Hz
fprintf('theory val freq %.0f�%.0f \n',to_theory.freq.val*1e-6,to_theory.freq.unc*1e-6)

to_old.wl.val=413.0938e-9;
to_old.wl.unc=sqrt(0.0009e-9^2+0.0020e-9^2);
[to_old.freq.val,to_old.freq.unc] = f2wl(to_old.wl.val,to_old.wl.unc); 

%% Generate data vectors
%to_val = [data.hwp.drift.to_val{1};data.qwp.drift.to_val{1};data.mix.drift.to_val{1}];%all our single scan tune-out values
%to_unc = [data.hwp.drift.to_val{2};data.qwp.drift.to_val{2};data.mix.drift.to_val{2}];%all corresponding uncertainties
to_val_for_polz=data.drift.to.val;
to_unc=data.drift.to.unc;
wlin=1./(to_unc.^2);

% to_val_for_polz=data.main.quad.to.val
% to_unc=data.main.quad.to.unc;
% wlin=1./(to_unc.^2);


%% polarisation model/data options
pol_opts.location = 'post';%post, pre_cen, pre_left, pre_right
pol_opts.predict_method = 'only_data';%'full_fit_pref_data','full_fit_pref_data','full_fit_only','only_data'; %obs (obsovation) fit (pertial fit) full_fit (fit with all parameters free)
                                        %'interp_only','interp_pref_data','interp_pref_interp'
                                        %'gauss_only','gauss_pref_data','gauss_pref_interp'
pol_opts.smoothing=3; %deg
pol_opts.wrap_hwp=0;

pol_opts.hwp=data.drift.wp.hwp;
pol_opts.qwp=data.drift.wp.qwp;
% 
% pol_opts.hwp=data.main.wp.hwp;
% pol_opts.qwp=data.main.wp.qwp;

pol_opts.add_noise=0; % used as a numerical way of propagating uncert
pol_model=pol_data_query(pol_opts);
pol_model_basic=pol_data_query_basic(pol_opts);

polz_theta_val=pol_model.theta.val;
polz_theta_ci=pol_model.theta.ci;
polz_v_val=pol_model.v.val;
polz_v_ci=pol_model.v.ci;
polz_cont_val=pol_model.cont.val;
polz_cont_ci=pol_model.cont.ci;

to_val_fit = to_val_for_polz;
wlin_fit = wlin./sum(wlin);
to_fit_val_offset=wmean(to_val_fit,wlin_fit); %should be very insensitive to this


% set up for the bootstrap wrapper
full_data_set=[polz_cont_val,polz_v_val,polz_theta_val,to_val_for_polz,to_unc];
full_data_set=num2cell(full_data_set,2);
% check that function gives the same answer as the above procedure
[wrap_fun_ans,detailed_out]=two_stage_two_dim_to_fit(full_data_set,to_fit_val_offset,1);
wrap_fun_ans=wrap_fun_ans*1e6+to_fit_val_offset;

to_scalar_minus_half_tensor=detailed_out.tsmht;
fit_mdl_full=detailed_out.fit;
fit_param_vals_full=fit_mdl_full.Coefficients.Estimate;
fit_param_vals_unc=fit_mdl_full.Coefficients.SE;

fprintf('to val model predict ...%.1f\n',(to_scalar_minus_half_tensor.val_predict*1e-6-725730000))
fprintf('to val model vals    ...%.1f\n',(to_scalar_minus_half_tensor.val*1e-6-725730000))
fprintf('to val wav    %.6f\n',f2wl(to_scalar_minus_half_tensor.val)*1e9)
fprintf('to unc model predict %.1f\n',(to_scalar_minus_half_tensor.unc_predict*1e-6))
fprintf('to unc model vals    %.1f\n',(to_scalar_minus_half_tensor.unc*1e-6))
fprintf('theta_k %.1f\n',mod(fit_param_vals_full(5),pi/2))
fprintf('angle offset %.3f\n',mod(fit_param_vals_full(4),pi))
fprintf('fit rmse %f\n',fit_mdl_full.RMSE)


%check that the tune out value is recovered with the correct query to the model function
fit_mdl_err=full_tune_out_polz_model(fit_param_vals_full,[-1,0,-fit_param_vals_full(4)])*1e6+to_fit_val_offset-to_scalar_minus_half_tensor.val;
if fit_mdl_err>eps
    error('cannot replicate prediction from wrapper function')
end

%%

% Do the bootstrap

% %detailed bootstrap
boot=bootstrap_se(@two_stage_two_dim_to_fit,full_data_set,...
    'opp_arguments',{to_fit_val_offset,1},...
    'plots',true,...
    'replace',true,...
    'plot_fig_name','TO boot MHz',...
    'samp_frac_lims',[0.10,1],...%[0.005,0.9]
    'num_samp_frac',3,... %20
    'num_samp_rep',3,... %1e2
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
fprintf('TO freq             %.1f�(%.0f(fit vals),%.0f(fit predict),%.0f�%.0f(boot)) MHz\n',...
    to_scalar_minus_half_tensor.val*1e-6,...
    to_scalar_minus_half_tensor.unc*1e-6,...
    to_scalar_minus_half_tensor.unc_predict*1e-6,...
    to_scalar_minus_half_tensor.unc_boot*1e-6,...
    to_scalar_minus_half_tensor.unc_unc_boot*1e-6)
fprintf('diff from TOV1      %.0f�%.0f MHz \n',(to_scalar_minus_half_tensor.val-to_old.freq.val)*1e-6,sqrt(to_scalar_minus_half_tensor.unc^2+to_old.freq.unc^2)*1e-6)
fprintf('diff from Theory    %.0f�%.0f MHz \n',...
    (to_scalar_minus_half_tensor.val-to_theory.freq.val)*1e-6,...
    sqrt(to_scalar_minus_half_tensor.unc^2+to_theory.freq.unc^2)*1e-6)
fprintf('TO wavelength       %.6f�%f nm \n',to_wav_val*1e9,to_wav_unc*1e9)

%% plot options

font_name='cmr10'; % ;
font_size_global=10;
font_size_label=18;

%% Find the polarization q paramter in the atomic frame
%we will compute the second Stokes parameter for the data
[polz_q_val,polz_q_ci]=Q_fun_unc(polz_cont_val,polz_theta_val,fit_param_vals_full(4),...
                                polz_cont_ci,polz_theta_ci,fit_param_vals_unc(4));

%%
% Full 2D plot of the fit in 2d Stokes space, with each scan shown
% we will plot in the second (Q) and fourth (v) Stokes parameters, 
% we rotate the measured Stokes parameters into the convinent basis found from the above fit



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
predict_predictor_input=cat(2,reshape(surf_polz_q,[],1),reshape(surf_polz_v,[],1),repmat(-fit_param_vals_full(4),numel(surf_polz_v),1));
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


scatter3(polz_v_val,polz_q_val,(to_val_for_polz-to_scalar_minus_half_tensor.val)*1e-6,'MarkerFaceColor',[0 .75 .75])
grid on
title('Full Stokes space representation')
xlabel('Fourth Stokes Parameter, V')
ylabel('Second Stokes Parameter, Q')
zlabel( sprintf('Tune out value-%.1f�%.1f (MHz)',to_scalar_minus_half_tensor.val*1e-6,to_scalar_minus_half_tensor.unc_boot*1e-6))
set(gcf,'color','w')
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1600, 900])
zlim([-8,10]*1e3)
view(-25,45)

stfig('2d TO residuals');
to_mdl_val=surf_mdl_fun(fit_mdl_full.Coefficients.Estimate,[polz_q_val,polz_v_val]);
to_mdl_val=to_mdl_val*1e6+to_fit_val_offset;
to_fit_resid=to_val_for_polz-to_mdl_val;
scatter3(polz_v_val,polz_q_val,to_fit_resid*1e-6,'MarkerFaceColor',[0 .75 .75])
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

%% Full 2D plot of the fit in 2d Stokes space, with each run binned into a single point
% we will plot in the second (Q) and fourth (v) Stokes parameters


colors_main=[[53,126,220];[33,188,44];[0,0,0]]./255;
lch=colorspace('RGB->LCH',colors_main(:,:));
lch(:,1)=lch(:,1)+20;
colors_detail=colorspace('LCH->RGB',lch);
%would prefer to use srgb_2_Jab here
color_shaded=colorspace('RGB->LCH',colors_main(3,:));
color_shaded(1)=50;
color_shaded=colorspace('LCH->RGB',color_shaded);


stfig('2d TO fit, binned');
scale_factor_freq=1e-9;
clf
surf_mdl=surface(surf_polz_v,surf_polz_q,(surf_samp_to-to_scalar_minus_half_tensor.val)*scale_factor_freq);
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
polz_vq=[polz_v_val,polz_q_val];
unique_polz_vq = unique(polz_vq, 'rows');
unique_polz_vq=unique_polz_vq(sum(isnan(unique_polz_vq),2)==0,:);


% use col_row_fun_mat to run the operation on rows of the matrix unique_polz_vq
% i select the elements of to_val_for_polz that corespond to matching the polz_vq matrix rows and thebinned_to_vq.to.valn
% compute the unc_wmean which returns a vector of weighted mean,weighted unc & unweighted std
iimax=size(unique_polz_vq,1);
binned_to_vq=[];
binned_to_vq.to=[];
binned_to_vq.to.val= nan(iimax,1);
binned_to_vq.to.ste= nan(iimax,1);
binned_to_vq.to.std= nan(iimax,1);
binned_to_vq.v=[];
binned_to_vq.v.val= nan(iimax,1);
binned_to_vq.v.ci= nan(iimax,2);
binned_to_vq.q=[];
binned_to_vq.q.val= nan(iimax,1);
binned_to_vq.q.ci= nan(iimax,2);

nan(iimax,9); % v,q,to,to_ste,to_std,v_ci+,v_ci-,q_unc
for ii=1:iimax
    target_vq=unique_polz_vq(ii,:);
    binned_to_vq.v.val(ii)=target_vq(1);
    binned_to_vq.q.val(ii)=target_vq(2);
    mask=all(polz_vq==target_vq,2);
    %polz_vq(mask,:) % dummy check 
    tmp_to_val_ste_std=unc_wmean_vec(to_val_for_polz(mask),to_unc(mask));
    binned_to_vq.to.val(ii)= tmp_to_val_ste_std(1);
    binned_to_vq.to.ste(ii)= tmp_to_val_ste_std(2);
    binned_to_vq.to.std(ii)= tmp_to_val_ste_std(3);
    binned_to_vq.v.ci(ii,:)=mean(polz_v_ci(mask,:),1);
    binned_to_vq.q.ci(ii,:)=mean(polz_q_ci(mask,:),1);
end



      
%scale
binned_to_vq_scaled=binned_to_vq;
binned_to_vq_scaled.to.val= (binned_to_vq_scaled.to.val-to_scalar_minus_half_tensor.val)*scale_factor_freq;
binned_to_vq_scaled.to.ste= binned_to_vq_scaled.to.ste*scale_factor_freq;
binned_to_vq_scaled.to.std= binned_to_vq_scaled.to.std*scale_factor_freq;

plot3d_ci_bars(binned_to_vq_scaled.v.val, binned_to_vq_scaled.q.val, binned_to_vq_scaled.to.val, ...
    binned_to_vq_scaled.v.ci(:,1),binned_to_vq_scaled.v.ci(:,2), ...
    binned_to_vq_scaled.q.ci(:,1),binned_to_vq_scaled.q.ci(:,2), ...
    -binned_to_vq_scaled.to.ste,binned_to_vq_scaled.to.ste, ...
    colors_main(1,:))
%plot3d_errorbars(polz_vq_to_val_unc_scaled(:,1), polz_vq_to_val_unc_scaled(:,2), polz_vq_to_val_unc_scaled(:,3), [], [], polz_vq_to_val_unc_scaled(:,4),colors_main(1,:));

scatter3(binned_to_vq_scaled.v.val, binned_to_vq_scaled.q.val, binned_to_vq_scaled.to.val,50,'o',...
   'MarkerEdgeColor',colors_main(1,:),'MarkerFaceColor',colors_detail(1,:),'LineWidth',2)
scatter3(0,-1,0,150,'D','MarkerFaceColor','k','MarkerEdgeColor','k')
grid on
%title('Full Stokes space representation')

set(gca,'FontSize',font_size_global,'FontName',font_name)
%xlabel('$\mathcal{V}$ ($4^{\mathrm{th}}$ Stokes parameter)','interpreter','latex','FontSize',font_size_label)
xlabel('$\mathcal{V}$','interpreter','latex','FontSize',font_size_label)
ylabel('$\mathcal{Q_{A}}$','interpreter','latex','FontSize',font_size_label) % ($2^{\mathrm{nd}}$ Stokes parameter)
zlabel('$\omega_{TO}(\mathcal{Q_{A}},\mathcal{V}) -\omega^{SMHT}_{TO}$ (GHz)','interpreter','latex','FontSize',font_size_label-2)
set(gcf,'Units','Pixels')
%set(gcf,'Position',[1068,355,676,453])
set(gca,'linewidth', 1.5)
set(gca,'TickLength',[0.02,0])
xticks([-1,-0.5,0,0.5,1])
yticks([-1,-0.5,0,0.5,1])
zticks(-8:4:10)
set(gcf,'color','w')
set(gcf,'Position',[1068,355,600,500])
view(-05-2*40,35)
%
%set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1600, 900])
zlim([-9,9])
drawnow
axis vis3d
%%
%view(-25,45)
for ii=0:2
    
    view(-10-ii*30,35)
    drawnow
    pause(0.1)
    export_fig(sprintf('.\\results\\20191112\\to_vq_dependence_%u.png',ii),'-a4','-r200')
end
hold off
axis vis3d

%%


%%
binned_to_vq_mdl_subtract=binned_to_vq;



stfig('2d TO residuals, binned');
to_mdl_val=surf_mdl_fun(fit_mdl_full.Coefficients.Estimate,[binned_to_vq_mdl_subtract.q.val,binned_to_vq_mdl_subtract.v.val]);
to_mdl_val=(to_mdl_val*1e6+to_fit_val_offset);

binned_to_vq_mdl_subtract.to.val=binned_to_vq_mdl_subtract.to.val-to_mdl_val;



scatter3(binned_to_vq_mdl_subtract.v.val, binned_to_vq_mdl_subtract.q.val, binned_to_vq_mdl_subtract.to.val,'MarkerFaceColor',[0 .75 .75])
hold on
plot3d_ci_bars(binned_to_vq_mdl_subtract.v.val, binned_to_vq_mdl_subtract.q.val, binned_to_vq_mdl_subtract.to.val, ...
    binned_to_vq_mdl_subtract.v.ci(:,1),binned_to_vq_mdl_subtract.v.ci(:,2),...
     [],[], ...
     -binned_to_vq_mdl_subtract.to.ste,binned_to_vq_mdl_subtract.to.ste,...
     'k');
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
%zlim([-1,1]*1e3)
view(-25,45)
xlabel('Fourth Stokes Parameter, V')
ylabel('Second Stokes Parameter, Q')
zlabel('Measurments-model (MHz)')



%% Plots of tune-out dependance on the individual Stokes parameters

% plot with Q, Stokes 2nd
samp_q=col_vec(linspace(-1.1,1.1,1e4));
predict_predictor_input=cat(2,samp_q,samp_q*0,repmat(-fit_param_vals_full(4),numel(samp_q),1));
[fit_mdl_val,fit_mdl_ci_curve]=predict(fit_mdl_full,predict_predictor_input,'Alpha',1-erf(1/sqrt(2)),'Prediction','Curve');
%'Prediction','Curve','Prediction','observation'
fit_mdl_val=fit_mdl_val*1e6+to_fit_val_offset;

fit_mdl_ci_curve=fit_mdl_ci_curve*1e6+to_fit_val_offset;
stfig('Q dep');
clf
patch([samp_q', fliplr(samp_q')],...
    ([fit_mdl_ci_curve(:,1)', fliplr(fit_mdl_ci_curve(:,2)')]-to_scalar_minus_half_tensor.val).*1e-6,...
    [1,1,1].*0.7,'EdgeColor','none')
hold on
plot(samp_q,(fit_mdl_val-to_scalar_minus_half_tensor.val)*1e-6,'k')
xlabel('Q value')
ylabel('Tune out -TOSMHT (MHz)')
title('V=0 extrapolation')
% now I want to plot every scan with its error bar that has been corrected onto V=0

shift_vq=v_correcting_shift(polz_v_val,polz_q_val,fit_mdl_full);

shifted_scaled_to_vals=(to_val_for_polz-to_scalar_minus_half_tensor.val+shift_vq(:,1)*1e6)*1e-6;

errorbar(polz_q_val,shifted_scaled_to_vals,to_unc*1e-6...
     ,'o','CapSize',0,'MarkerSize',5,'Color','r',...
     'LineWidth',1.5);
%ylim([min(shifted_scaled_to_vals),max(shifted_scaled_to_vals)])
ylim([-0.5,0.5]*1e4)
xlim([-1.1,1.1])

% do the same for V, Stokes 4th
samp_v=col_vec(linspace(-1.1,1.1,1e4));
predict_predictor_input=cat(2,samp_v*0-1,samp_v,repmat(-fit_param_vals_full(4),numel(samp_v),1));
[fit_mdl_val,fit_mdl_ci_curve]=predict(fit_mdl_full,predict_predictor_input,'Alpha',1-erf(1/sqrt(2)),'Prediction','Curve');
%'Prediction','Curve','Prediction','observation'
fit_mdl_val=fit_mdl_val*1e6+to_fit_val_offset;
fit_mdl_ci_curve=fit_mdl_ci_curve*1e6+to_fit_val_offset;
stfig('V dep');
clf
patch([samp_v', fliplr(samp_v')],...
    ([fit_mdl_ci_curve(:,1)', fliplr(fit_mdl_ci_curve(:,2)')]-to_scalar_minus_half_tensor.val).*1e-6,...
    [1,1,1].*0.7,'EdgeColor','none')
hold on
plot(samp_v,(fit_mdl_val-to_scalar_minus_half_tensor.val)*1e-6,'k')
xlabel('V value')
ylabel('Tune out -TOSMHT(MHz)')
title('Q=-1 extrapolation')
shifted_scaled_to_vals=(to_val_for_polz-to_scalar_minus_half_tensor.val+shift_vq(:,2)*1e6)*1e-6;

errorbar(polz_v_val,shifted_scaled_to_vals,to_unc*1e-6...
     ,'o','CapSize',0,'MarkerSize',5,'Color','r',...
     'LineWidth',1.5);
%ylim([min(shifted_scaled_to_vals),max(shifted_scaled_to_vals)]) 
ylim([-1,0.8].*1e4) 
xlim([-1.1,1.1])
 

%% Plot the V,Q plots with binned data

font_size_global=18;
font_size_label=19;



colors_main=[[53,126,220];[33,188,44];[0,0,0]]./255;
lch=colorspace('RGB->LCH',colors_main(:,:));
lch(:,1)=lch(:,1)+20;
colors_detail=colorspace('LCH->RGB',lch);
%would prefer to use srgb_2_Jab here
color_shaded=colorspace('RGB->LCH',colors_main(3,:));
color_shaded(1)=50;
color_shaded=colorspace('LCH->RGB',color_shaded);


% plot with Q, Stokes 2nd
scale_factor_freq=1e-9;
samp_q=col_vec(linspace(-1.1,1.1,1e4));
predict_predictor_input=cat(2,samp_q,samp_q*0,repmat(-fit_param_vals_full(4),numel(samp_q),1));
[fit_mdl_val,fit_mdl_ci_curve]=predict(fit_mdl_full,predict_predictor_input,'Alpha',1-erf(1/sqrt(2)),'Prediction','Curve');
%'Prediction','Curve','Prediction','observation'
fit_mdl_val=fit_mdl_val*1e6+to_fit_val_offset;
fit_mdl_ci_curve=fit_mdl_ci_curve*1e6+to_fit_val_offset;

[~,fit_mdl_ci_obs]=predict(fit_mdl_full,predict_predictor_input,'Alpha',1-erf(1/sqrt(2)),'Prediction','observation');
fit_mdl_ci_obs=fit_mdl_ci_obs*1e6+to_fit_val_offset;
stfig('Q dep,binned');
clf
patch([samp_q', fliplr(samp_q')],...
    ([fit_mdl_ci_curve(:,1)', fliplr(fit_mdl_ci_curve(:,2)')]-to_scalar_minus_half_tensor.val).*scale_factor_freq,...
    [1,1,1].*0.7,'EdgeColor','none')

patch([samp_q', fliplr(samp_q')],...
    ([fit_mdl_ci_obs(:,1)', fliplr(fit_mdl_ci_obs(:,2)')]-to_scalar_minus_half_tensor.val).*scale_factor_freq,...
    [1,1,1].*0.5,'EdgeColor','none')

hold on
plot(samp_q,(fit_mdl_val-to_scalar_minus_half_tensor.val)*scale_factor_freq,'k')

%
% now I want to plot every scan with its error bar that has been corrected onto V=0
% caluate the frequnecy shift for each binv,q
shift_vq=v_correcting_shift(binned_to_vq.v.val,binned_to_vq.q.val,fit_mdl_full);
shifted_scaled_to_vals=(binned_to_vq.to.val-to_scalar_minus_half_tensor.val+shift_vq(:,1)*1e6)*scale_factor_freq;

deriv_to_with_v=derivest(@(x) predict(fit_mdl_full,[0,x,-fit_param_vals_full(4)]),1,'Vectorized','no');
deriv_to_with_v=deriv_to_with_v*1e6;

combined_to_unc=(binned_to_vq.v.ci*deriv_to_with_v).^2+([-binned_to_vq.to.ste,binned_to_vq.to.ste]).^2;
combined_to_unc=[-sqrt(combined_to_unc(:,1)),sqrt(combined_to_unc(:,2))];
%

% version thin bars for polz
errorbar(binned_to_vq.q.val,shifted_scaled_to_vals, ...
    combined_to_unc(:,1)*scale_factor_freq,combined_to_unc(:,1)*scale_factor_freq,...
    binned_to_vq.q.ci(:,1),binned_to_vq.q.ci(:,1),...
     'o','CapSize',0,'Marker','none','Color',colors_detail(1,:),...
     'LineWidth',1.5);
errorbar(binned_to_vq.q.val,shifted_scaled_to_vals, ...
        binned_to_vq.to.ste*scale_factor_freq,binned_to_vq.to.ste*scale_factor_freq,...
    'o','CapSize',0,'MarkerSize',5,'Color',colors_main(1,:),...
    'MarkerFaceColor',colors_detail(1,:),'LineWidth',2.8);


plot(-1,0,'xr')
%ylim([min(shifted_scaled_to_vals),max(shifted_scaled_to_vals)])
ylim([-2.5,1.2]*1)
xlim([-1.1,1.1])
hold off
box on
set(gca,'FontSize',font_size_global,'FontName',font_name)
xlabel('$\mathcal{Q_{A}}$ ($2^{\mathrm{nd}}$ Stokes parameter)','interpreter','latex','FontSize',font_size_label)
ylabel('$f_{TO}(\mathcal{Q_{A}},\mathcal{V}=0) -f_{TO}(-1,0)$ (GHz)','interpreter','latex','FontSize',font_size_label)
set(gcf,'Units','Pixels')
set(gcf,'Position',[1068,355,676,453])
set(gca,'linewidth', 1.5)
set(gca,'TickLength',[0.02,0])
xticks([-1,-0.5,0,0.5,1])
yticks(-2:1:1)

%%

stfig('Q dep,binned,residuals');
subplot(2,1,1)
samp_q=binned_to_vq.q.val;
predict_predictor_input=cat(2,samp_q,samp_q*0,repmat(-fit_param_vals_full(4),numel(samp_q),1));
fit_mdl_val=predict(fit_mdl_full,predict_predictor_input,'Alpha',1-erf(1/sqrt(2)),'Prediction','Curve');
%'Prediction','Curve','Prediction','observation'
fit_mdl_val=fit_mdl_val*1e6+to_fit_val_offset;
shifted_scaled_to_vals=(binned_to_vq.to.val+shift_vq(:,1)*1e6-fit_mdl_val)*1e-6;
errorbar(samp_q,shifted_scaled_to_vals,...
    combined_to_unc(:,1)*1e-6,combined_to_unc(:,2)*1e-6...
     ,'o','CapSize',0,'Marker','none','Color','g',...
     'LineWidth',1.5);
hold on
% errorbar(samp_q,shifted_scaled_to_vals,polz_vq_to_val_unc(:,4)*1e-6...
%      ,'o','CapSize',0,'MarkerSize',5,'Color','r',...
%      'LineWidth',1.5);
% hold off
subplot(2,1,2)
% now find this in terms of standard deviations
% need to select the right CI depending if pos or neg
resid_neg_mask=shifted_scaled_to_vals<=0;
% get the ci corresponding to the sign of the resid
signed_ci=combined_to_unc(:,1)*nan;
signed_ci(resid_neg_mask)=combined_to_unc(resid_neg_mask,1);
signed_ci(~resid_neg_mask)=combined_to_unc(~resid_neg_mask,2);
residuals_num_ste=shifted_scaled_to_vals./(signed_ci*1e-6);
fprintf('standard deviation of errors/ci %f \n',std(residuals_num_ste))
plot(samp_q,residuals_num_ste,'ok')
fprintf('sd of (residuals/ste) in Q dep %.1f \n',std(residuals_num_ste))
xlabel('$\mathcal{Q_{A}}$, ($2^{\mathrm{nd}}$ Stokes parameter )') %font_size_label
ylabel('Error From Model (MHz)')
set(gca,'FontSize',font_size_global,'FontName',font_name)





%%

% do the same for V, Stokes 4th
% extrapolate to Q=0
samp_v=col_vec(linspace(-1.1,1.1,1e4));
scale_factor_freq=1e-9;
predict_predictor_input=cat(2,samp_v*0,samp_v,repmat(-fit_param_vals_full(4),numel(samp_v),1));
[fit_mdl_val,fit_mdl_ci_curve]=predict(fit_mdl_full,predict_predictor_input,'Alpha',1-erf(1/sqrt(2)),'Prediction','Curve');
%'Prediction','Curve','Prediction','observation'
fit_mdl_val=fit_mdl_val*1e6+to_fit_val_offset;
fit_mdl_ci_curve=fit_mdl_ci_curve*1e6+to_fit_val_offset;
stfig('V dep,binned');
clf
patch([samp_v', fliplr(samp_v')],...
    ([fit_mdl_ci_curve(:,1)', fliplr(fit_mdl_ci_curve(:,2)')]-to_scalar_minus_half_tensor.val).*scale_factor_freq,...
    [1,1,1].*0.7,'EdgeColor','none')
hold on
plot(samp_v,(fit_mdl_val-to_scalar_minus_half_tensor.val)*scale_factor_freq,'k')
%xlabel('$4^{th}$ Stokes parameter,\bfv\rm','interpreter','latex')
xlabel('$\mathcal{V}$ ($4^{\mathrm{th}}$ Stokes parameter)','interpreter','latex')
ylabel('$\omega_{TO}(\mathcal{Q_{A}}=0,\mathcal{V}) -\omega^{SMHT}_{TO}$ (GHz)','interpreter','latex')
%title('Q=-1 extrapolation')

shifted_scaled_to_vals=(polz_vq_to_val_unc(:,3)-to_scalar_minus_half_tensor.val+shift_vq(:,2)*1e6)*scale_factor_freq;
errorbar(polz_vq_to_val_unc(:,1),shifted_scaled_to_vals,polz_vq_to_val_unc(:,5)*scale_factor_freq...
     ,'o','CapSize',0,'Marker','none','Color',colors_detail(1,:),...
     'LineWidth',1.5);
errorbar(polz_vq_to_val_unc(:,1),shifted_scaled_to_vals,polz_vq_to_val_unc(:,4)*scale_factor_freq,...
    'o','CapSize',0,'MarkerSize',5,'Color',colors_main(1,:),...
    'MarkerFaceColor',colors_detail(1,:),'LineWidth',2.8);
hold off
%ylim([min(shifted_scaled_to_vals),max(shifted_scaled_to_vals)]) 
ylim([-8.5,8.5].*1) 
xlim([-1.1,1.1])
xticks([-1,-0.5,0,0.5,1])
yticks(-8:4:8)
box on
set(gca,'FontSize',font_size_global,'FontName',font_name)
set(gca,'linewidth', 1.5)
set(gca,'TickLength',[0.02,0])
set(gcf,'Units','Pixels')
%set(gcf,'Position',[1068,355,676,453])




stfig('V dep,binned,residuals');
subplot(2,1,1)
samp_v=polz_vq_to_val_unc(:,1);
predict_predictor_input=cat(2,samp_v*0-1,samp_v,repmat(-fit_param_vals_full(4),numel(samp_v),1));
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


function [polz_q_val,polz_q_unc]=Q_fun_unc(cont_val,theta_val,phi_val,...
                                cont_ci,theta_ci,phi_unc)
        polz_q_val=Q_fun(cont_val,theta_val,phi_val);
        if nargout>1
            % propagate uncertianty as tolerances
            % this isnt a good way to do this but it is quick to implement
            % i can check this with numerical experiments later
            % want all possible combinations of 3 values drawn from [-1,0,1] with replacements
            signs=unique(nchoosek(repelem([-1,0,1],3), 3),'rows'); % this is a inefficeint way to get this, but all i can find with inbuilts
            ans_range=zeros(size(polz_q_val,1),size(signs,1));
            for ii=1:size(signs,1)
                if signs(ii,1)==1
                    cont_delt=cont_ci(:,2);
                elseif signs(ii,1)==-1
                    cont_delt=cont_ci(:,1);
                else 
                    cont_delt=0;
                end
                
                if signs(ii,2)==1
                    theta_delt=theta_ci(:,2);
                elseif signs(ii,2)==-1
                    theta_delt=theta_ci(:,1);
                else 
                    theta_delt=0;
                end
                
                phi_delt=signs(ii,3)*phi_unc;

                ans_range(:,ii)=Q_fun(cont_val+cont_delt,theta_val+theta_delt,phi_val+phi_delt);
            end
            ans_ci=[nanmin(ans_range,[],2),nanmax(ans_range,[],2)];
            ans_ci=ans_ci-repmat(polz_q_val,[1,2]);
            %[ans_ci,polz_q_val]
            polz_q_unc=ans_ci;
        end
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
    
    predict_predictor_input=cat(2,q_in*0,v_in,repmat(-mdl.Coefficients.Estimate(4),numel(v_in),1));
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


% find the tune out value two ways the first from the fit model asking what the freq will be at the appropriate
% polarization parameters
[tsmht_details.val_predict,tsmht_details.unc_predict]=...
    predict(fit_mdl_full,[-1,0,-fit_vals_full(4)],'Alpha',1-erf(1/sqrt(2)),'Prediction','Curve');
tsmht_details.val_predict=tsmht_details.val_predict*1e6+fit_offset;
tsmht_details.unc_predict=1/2*range(tsmht_details.unc_predict)*1e6;

% then the next method uses the fit parameters and calculates the freq using those fit parameters
% the main reason to do this two ways is to look at the uncert

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

%V = zeros(numel(to_val),1);%fourth Stokes parameter
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
% polz_q = Q_fun(polz_cont,polz_theta,fit_vals_full(4));%second Stokes parameter, relative to the vector that is otrhogonal to the B feild and beam axis
% % TO CHECK
% %Q_run = Q_fun(polz_cont,polz_theta,fit_vals(4));%second Stokes parameter for each run
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
% disp_config.fig_name = 'TO fit Stokes 4';
% disp_config.plot_title='';
% disp_config.x_label = 'Fourth Stokes Parameter, V';
% disp_config.plot_offset.val=fit_vals(1)./1e6;
% disp_config.plot_offset.unc=fit_uncs(1)./1e6;
% %plot_binned_nice(disp_config,x_dat,y_dat,weights,fit_mdl)
% plot_binned_nice(disp_config,polz_v,tens_corr_to./1e6,wlin,fit_mdl_v)
% 
% disp_config.fig_name = 'TO fit Stokes 2';
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