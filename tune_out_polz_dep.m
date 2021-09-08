clear all


%% Exp parameters
aom_freq =0; %aom offset in Hz, as of 20190528 offset is calculated in main_trap_freq

%%
% addpath('../Core_BEC_Analysis/') %add the path to set_up_project_path
addpath('./data/') %add the path to set_up_project_path
% set_up_project_path


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
root_data_dir='E:\Data\Tuneout\to_main_data';
files = dir(root_data_dir);
files=files(3:end);
% Get a logical vector that tells which is a directory.
dir_mask = [files.isdir];
folders=files(dir_mask);
folders=arrayfun(@(x) fullfile(root_data_dir,x.name),folders,'UniformOutput' ,false);
folders;
loop_config.dir=folders;
%%

data = load_pocessed_to_data(loop_config);

%% dist of atom number
qe=0.09;
[mean_atomnum,std_atomnum]=pooled_mean_and_std(data.drift.atom_num_probe.val,data.drift.atom_num_probe.std,data.drift.probe_shots);
mean_atomnum=mean_atomnum/qe;
std_atomnum=std_atomnum/qe;

smooth_hist(data.drift.atom_num_probe.val/qe,'sigma',1.5e4);
hold on
smooth_hist(data.drift.atom_num_probe.val/qe,'sigma',3e3);
hold off

%%
%save('./data/20191119_imported_data_for_tune_out_polz_dep.mat')
%load('./data/20190611_imported_data_for_tune_out_polz_dep.mat')

%% Theory
% TO val from https://journals.aps.org/pra/abstract/10.1103/PhysRevA.93.052516, 413.0859(4)
% half tensor shift from \theta_p=54.74 \lambda_TO-> 413.083876e-9
% half tensor shift from \theta_p=90 \lambda_TO-> 413.082896e-9
to_theory=[];
to_theory.freq.val=725736430e6;
to_theory.freq.unc=50;
[to_theory.wl.val,to_theory.wl.unc] = f2wl(to_theory.freq.val,to_theory.freq.unc); %best theory guess in Hz
fprintf('theory val freq %.0f±%.0f \n',to_theory.freq.val*1e-6,to_theory.freq.unc*1e-6)

to_old.wl.val=413.0938e-9;
to_old.wl.unc=sqrt(0.0009e-9^2+0.0020e-9^2);
[to_old.freq.val,to_old.freq.unc] = f2wl(to_old.wl.val,to_old.wl.unc); 

% % Generate data vectors
%to_val = [data.hwp.drift.to_val{1};data.qwp.drift.to_val{1};data.mix.drift.to_val{1}];%all our single scan tune-out values
%to_unc = [data.hwp.drift.to_val{2};data.qwp.drift.to_val{2};data.mix.drift.to_val{2}];%all corresponding uncertainties
to_val_for_polz=data.drift.to.val;
to_unc=data.drift.to.unc;
wlin=1./(to_unc.^2);

%% polarisation model/data options
pol_opts.location = 'post';%post, pre_cen, pre_left, pre_right
pol_opts.predict_method = 'only_data';%'full_fit_pref_fit','full_fit_pref_data','full_fit_only','only_data'; %obs (obsovation) fit (pertial fit) full_fit (fit with all parameters free)
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
to_fit_val_offset=wmean(to_val_fit,wlin_fit); %should be very insensitive to this


% set up for the bootstrap wrapper
full_data_set=[polz_cont,polz_v,polz_theta,to_val_for_polz,to_unc];
full_data_set=num2cell(full_data_set,2);
% check that function gives the same answer as the above procedure
[wrap_fun_ans,detailed_out]=two_stage_two_dim_to_fit(full_data_set,to_fit_val_offset,1);
wrap_fun_ans=wrap_fun_ans*1e6+to_fit_val_offset;

% % HMMM.
% Well - what if we looked at the delta (pred) as a fn of parameters (which
% we should see in Jacobian...? nah; that's change in
% predicted-out-X-values wrt each parameter


to_scalar_minus_half_tensor=detailed_out.tsmht;
fit_mdl_full=detailed_out.fit;
fit_vals_full=fit_mdl_full.Coefficients.Estimate;

fprintf('to val model predict ...%.1f\n',(to_scalar_minus_half_tensor.val_predict*1e-6-725730000))
fprintf('to val model vals    ...%.1f\n',(to_scalar_minus_half_tensor.val*1e-6-725730000))
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
%%



%%

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

font_name='cmr10'; % ;
font_size_global=10;
font_size_label=18;


%%
% Full 2D plot of the fit in 2d Stokes space, with each scan shown
% we will plot in the second (Q) and fourth (v) Stokes parameters, 
% we rotate the measured Stokes parameters into the convinent basis found from the above fit

%we will compute the second Stokes parameter for the data
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
title('Full Stokes space representation')
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
polz_vq_to_val_unc_scaled(:,3)=(polz_vq_to_val_unc(:,3)-to_scalar_minus_half_tensor.val)*scale_factor_freq;
polz_vq_to_val_unc_scaled(:,4)=polz_vq_to_val_unc(:,4)*scale_factor_freq; %scale the unc
polz_vq_to_val_unc_scaled(:,5)=polz_vq_to_val_unc(:,5)*scale_factor_freq; %scale the std

plot3d_errorbars(polz_vq_to_val_unc_scaled(:,1), polz_vq_to_val_unc_scaled(:,2), polz_vq_to_val_unc_scaled(:,3), [], [], polz_vq_to_val_unc_scaled(:,4),colors_main(1,:));

scatter3(polz_vq_to_val_unc_scaled(:,1), polz_vq_to_val_unc_scaled(:,2), polz_vq_to_val_unc_scaled(:,3),50,'o',...
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
% for ii=0:2
%     
%     view(-10-ii*30,35)
%     drawnow
%     pause(0.1)
%     export_fig(sprintf('.\\results\\20191112\\to_vq_dependence_%u.png',ii),'-a4','-r200')
% end
% hold off
% axis vis3d

%%


%%
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



%% Plots of tune-out dependance on the individual Stokes parameters

% plot with Q, Stokes 2nd
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

% do the same for V, Stokes 4th
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
samp_q=col_vec(linspace(-1.1,1.1,1e2));
predict_predictor_input=cat(2,samp_q,samp_q*0,repmat(-fit_vals_full(4),numel(samp_q),1));
[fit_mdl_val,fit_mdl_ci]=predict(fit_mdl_full,predict_predictor_input,'Alpha',1-erf(1/sqrt(2)),'Prediction','Curve');
[obs_fit_mdl_val,obs_fit_mdl_ci]=predict(fit_mdl_full,predict_predictor_input,'Alpha',0.05,'Prediction','observation','Simultaneous','true');
%'Prediction','Curve','Prediction','observation'
fit_mdl_val=fit_mdl_val*1e6+to_fit_val_offset;
fit_mdl_ci=fit_mdl_ci*1e6+to_fit_val_offset;

obs_fit_mdl_val=obs_fit_mdl_val*1e6+to_fit_val_offset;
obs_fit_mdl_ci=obs_fit_mdl_ci*1e6+to_fit_val_offset;


stfig('Q dep,binned');
clf
% subplot(2,1,1)
hold on
patch([samp_q', fliplr(samp_q')],...
    ([obs_fit_mdl_ci(:,1)', fliplr(obs_fit_mdl_ci(:,2)')]-to_scalar_minus_half_tensor.val).*scale_factor_freq,...
    [1,0,1].*0.6,'EdgeColor','none','FaceAlpha',0.5)
%
% subplot(2,1,2)
hold on
patch([samp_q', fliplr(samp_q')],...
    ([fit_mdl_ci(:,1)', fliplr(fit_mdl_ci(:,2)')]-to_scalar_minus_half_tensor.val).*scale_factor_freq,...
    [0,0,1],'EdgeColor','none','FaceAlpha',0.2)



plot(samp_q,(fit_mdl_val-to_scalar_minus_half_tensor.val)*scale_factor_freq,'k')

%title('v=0 extrapolation')
% now I want to plot every scan with its error bar that has been corrected onto V=0

shift_vq=v_correcting_shift(polz_vq_to_val_unc(:,1),polz_vq_to_val_unc(:,2),fit_mdl_full);
shifted_scaled_to_vals=(polz_vq_to_val_unc(:,3)-to_scalar_minus_half_tensor.val+shift_vq(:,1)*1e6)*scale_factor_freq;

errorbar(polz_vq_to_val_unc(:,2),shifted_scaled_to_vals,polz_vq_to_val_unc(:,5)*scale_factor_freq...
     ,'o','CapSize',0,'Marker','none','Color',colors_detail(1,:),...
     'LineWidth',1.5);
errorbar(polz_vq_to_val_unc(:,2),shifted_scaled_to_vals,polz_vq_to_val_unc(:,4)*scale_factor_freq,...
    'o','CapSize',0,'MarkerSize',5,'Color',colors_main(1,:),...
    'MarkerFaceColor',colors_detail(1,:),'LineWidth',2.8);
ylim([min(shifted_scaled_to_vals),max(shifted_scaled_to_vals)])
ylim([-2.5,1.2]*1)
xlim([-1.1,1.1])
hold off
box on
set(gca,'FontSize',font_size_global,'FontName',font_name)
xlabel('$\mathcal{Q_{A}}$ ($2^{\mathrm{nd}}$ Stokes parameter)','interpreter','latex','FontSize',font_size_label)
% ylabel('$\omega_{TO}(\mathcal{Q_{A}},\mathcal{V}=0) -\omega^{SMHT}_{TO}$ (GHz)','interpreter','latex','FontSize',font_size_label)
set(gcf,'Units','Pixels')
%set(gcf,'Position',[1068,355,676,453])
set(gca,'linewidth', 1.5)
set(gca,'TickLength',[0.02,0])
xticks([-1,-0.5,0,0.5,1])
% yticks(-2:1:1)



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
xlabel('$\mathcal{Q_{A}}$, ($2^{\mathrm{nd}}$ Stokes parameter )') %font_size_label
ylabel('Error From Model (MHz)')
set(gca,'FontSize',font_size_global,'FontName',font_name)







% do the same for V, Stokes 4th
% extrapolate to Q=0
samp_v=col_vec(linspace(-1.1,1.1,1e4));
scale_factor_freq=1e-9;
predict_predictor_input=cat(2,samp_v*0,samp_v,repmat(-fit_vals_full(4),numel(samp_v),1));
[fit_mdl_val,fit_mdl_ci]=predict(fit_mdl_full,predict_predictor_input,'Alpha',1-erf(1/sqrt(2)),'Prediction','Curve');
[obs_fit_mdl_val,obs_fit_mdl_ci]=predict(fit_mdl_full,predict_predictor_input,'Alpha',.05,'Prediction','observation','simultaneous','true');
%'Prediction','Curve','Prediction','observation'
fit_mdl_val=fit_mdl_val*1e6+to_fit_val_offset;
fit_mdl_ci=fit_mdl_ci*1e6+to_fit_val_offset;
obs_fit_mdl_val=obs_fit_mdl_val*1e6+to_fit_val_offset;
obs_fit_mdl_ci=obs_fit_mdl_ci*1e6+to_fit_val_offset;
stfig('V dep,binned');
clf
hold on
patch([samp_v', fliplr(samp_v')],...
    ([fit_mdl_ci(:,1)', fliplr(fit_mdl_ci(:,2)')]-to_scalar_minus_half_tensor.val).*scale_factor_freq,...
    [1,1,1].*0.7,'EdgeColor','none','FaceAlpha',0.6)
patch([samp_v', fliplr(samp_v')],...
    ([obs_fit_mdl_ci(:,1)', fliplr(obs_fit_mdl_ci(:,2)')]-to_scalar_minus_half_tensor.val).*scale_factor_freq,...
    [0,0,1],'EdgeColor','none','FaceAlpha',0.6)

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

%% linear model?

% Issue with prior approach: Assumed linear in polz_q but did not account
% for transformation into atomic frame. This takes the form
% Q_A = fun(inputs,theta);
% where we cannot precisely determine theta_L. 
% Instead we can optimize theta_L with respect to the goodness of fit


% to try; compute the linear model coefs from the NL model. How close are
% they?

cli_header('Attempting to optimize theta_L...');

use_vals = ~any(isnan(polz_vq),2);
l_fit_offset = 0*to_theory.freq.val;
alpha_use = .05;
ci_use = 1 - alpha_use;
fit_scale_fac = 1;
to_val_for_lfit = (to_val_fit-l_fit_offset)/fit_scale_fac;




opt_theta_fit_call = @(x) opt_theta_eval(polz_cont,polz_theta,polz_v,to_val_for_lfit,to_unc,x);
theta_L = fminunc(opt_theta_fit_call,2)
Q_A=Q_fun(polz_cont,polz_theta,theta_L);
full_cov =fit_mdl_full.CoefficientCovariance;
[~,opt_out] = opt_theta_eval(polz_cont,polz_theta,polz_v,to_val_for_lfit,to_unc,theta_L);
cli_header(1,'Optimized value %.3f',theta_L);
lfit_init = opt_out.fit;
% ah but the fit assumes one also has the Q transformation right? so it
% just looks like a NLM. Or one could go
lfit_eval = @(polz_cont,theta_min,polz_v) lfit_init([Q_fun(polz_cont,theta_min,theta_L),polz_v]);
gof = opt_out.gof;
out = opt_out.out;
q_A_nl = Q_fun(polz_cont,polz_theta, fit_mdl_full.Coefficients.Estimate(4));
q_A = opt_out.predictors(:,1);
v_A = opt_out.predictors(:,2);
predictors = opt_out.predictors;
observations = opt_out.observations;


weights = to_unc(use_vals).^-2;
[lfit_fin,lfit_gof,lfit_out] = fit([q_A,v_A],to_val_for_lfit(use_vals),'poly11');
% %
pred_vals=lfit_fin([q_A,v_A]);
pred_PI = predint(lfit_fin,[q_A,v_A],erf(1/sqrt(2)),'observation','off')-pred_vals;
pred_std = mean(abs(pred_PI),2);
chi_square = sum((lfit_out.residuals./pred_std).^2);
chi_square_per_dof = sum((lfit_out.residuals./pred_std).^2)/lfit_gof.dfe;


% fds = cell2mat(full_data_set);
% chi_pred = [polz_q,polz_v,-fit_mdl_full.Coefficients.Estimate(4)*ones(size(polz_q))];
% nmask = ~isnan(polz_v);
% chi_pred = chi_pred(nmask,:);
% [pred_nl, pred_nl_ci] = predict(fit_mdl_full,chi_pred,'Alpha',1-erf(1/sqrt(2)),'Prediction','observation')
% pred_nl = pred_nl*1e6+to_fit_val_offset
% pred_nl_ci = pred_nl_ci*1e6+to_fit_val_offset-pred_nl
% pred_std = mean(abs(pred_nl_ci),2);
% nlm_residuals = pred_nl-to_val_for_polz(~isnan(fds(:,1)));
% chi_squared = sum(nlm_residuals.^2./pred_std.^2);
% dof = lfit_out.numobs - lfit_out.numparam;
% chi_squared_per_dof = chi_squared/dof
% %


% 1 sigma?

l_CI = confint(lfit_fin);
lfit_2 = fitlm(predictors,observations,'linear');

cli_header('Linear fit results');
lfit_cov = lfit_2.CoefficientCovariance;
% construct coefs from NL for comparison
% a_reduced_coef = fit_mdl_full.Coefficients.Estimate(1) - ...
%     fit_mdl_full.Coefficients.Estimate(3)*(0.5*0.75*(sin(fit_mdl_full.Coefficients.Estimate(5))).^2);
% b_reduced_coef = 0.5*fit_mdl_full.Coefficients.Estimate(2)*cos(fit_mdl_full.Coefficients.Estimate(5));
% c_reduced_coef = 0.75*fit_mdl_full.Coefficients.Estimate(3)*(sin(fit_mdl_full.Coefficients.Estimate(5))).^2;
% tranformed_coefs = [c_reduced_coef,b_reduced_coef,c_reduced_coef];




QV_pred = [-1,0];
pred_TO = lfit_fin(QV_pred)*fit_scale_fac + l_fit_offset;
pred_CI = predint(lfit_fin,QV_pred,ci_use,'functional','off')*fit_scale_fac+l_fit_offset;
pred_PI = predint(lfit_fin,QV_pred,ci_use,'observation','off')*fit_scale_fac+l_fit_offset;


Q_plot = linspace(-1.1,1.1,100);
V_plot = linspace(-1.1,1.1,100);
[Q,V] = meshgrid(Q_plot,V_plot);
pred_surf = fit_scale_fac*lfit_fin(Q,V)+l_fit_offset;

% will also want to look at prediction intervals. So for all the XY
% pairs...
Q_col=reshape(Q,numel(Q),1);
V_col=reshape(V,numel(V),1);
surf_PI=predint(lfit_fin,[Q_col,V_col])*1e6+l_fit_offset;


fprintf('\n');
cli_header('Linear fit done:');
cli_header(1,'%.1f pc CI for fit parameters:',100*ci_use);
cli_header(2,'%.f (%.f,%.f)',lfit_fin.p00,l_CI(:,1)')
cli_header(2,'%.f (%.f,%.f)',lfit_fin.p10,l_CI(:,2)')
cli_header(2,'%.f (%.f,%.f)',lfit_fin.p01,l_CI(:,3)')
cli_header(1,'Predicted f_TO(-1,0): %.f(%.f,%.f) MHz',1e-6*pred_TO,1e-6*(pred_CI-pred_TO));
cli_header(2,'Using non-simultaneous prediction CI' );
cli_header(1,'Prediction interval: (%.f,%.f) MHz',1e-6*(pred_PI-pred_TO));
cli_header(2,'Using non-simultaneous prediction CI' );
cli_header(1,'Differs from prior method prediction %.f (%.f) by %.f MHz',(to_scalar_minus_half_tensor.val)/1e6,(to_scalar_minus_half_tensor.unc_predict)/1e6,(to_scalar_minus_half_tensor.val-pred_TO)/1e6);


%well.. as a first step 

stfig('lfit tests');
clf
hold on
scatter3(predictors(:,1),predictors(:,2),scale_factor_freq*(observations-pred_TO),50,...
    'MarkerEdgeColor','k',...
        'MarkerFaceColor',[0 .75 .75])
% scatter3(-1,0,0,...
%     50,'MarkerEdgeColor','k',...
%         'MarkerFaceColor',[0.75 0 .75 ])
scatter3(polz_q(use_vals),predictors(:,2),scale_factor_freq*(observations-pred_TO),50,...
    'MarkerEdgeColor','k',...
        'MarkerFaceColor',[0 .75 0])
surf_mdl=surface(surf_polz_q,surf_polz_v,(surf_samp_to-to_scalar_minus_half_tensor.val)*scale_factor_freq);

% pi_l = surf(Q,V,scale_factor_freq*(reshape(surf_PI(:,1),size(Q))-pred_TO),...
%     'FaceColor',0.3*[1,1,1]);
% set(pi_l, 'FaceAlpha', 0.3)
% pi_u = surf(Q,V,scale_factor_freq*(reshape(surf_PI(:,2),size(Q))-pred_TO),...
%     'FaceColor',0.3*[1,1,1]);
% set(pi_u, 'FaceAlpha', 0.3)
% pi_l.CData = 10*ones(size(pred_surf));
% pi_u.CData = 10*ones(size(pred_surf));

surf_lin=surf(Q,V,scale_factor_freq*(pred_surf-pred_TO),'EdgeColor','none');

legend('$Q_A$ data','$Q_L$ data','NL fit','L fit')

set(surf_lin, 'FaceAlpha', 0.7)
surf_lin.CData = scale_factor_freq*(pred_surf-pred_TO);
% colormap(viridis)

ylabel('$\mathcal{V}$','interpreter','latex','FontSize',font_size_label)
xlabel('$\mathcal{Q_{A}}$','interpreter','latex','FontSize',font_size_label) % ($2^{\mathrm{nd}}$ Stokes parameter)
zlabel('$\omega_{TO}(\mathcal{Q_{A}},\mathcal{V}) -\omega^{SMHT}_{TO}$ (GHz)','interpreter','latex','FontSize',font_size_label-2)

shading interp
% %

% %


% last step: include the full nonlinear model in the interval analysis
% detailed_out has the goods



% How to project to 1D?
% Some options;
%     As above, predict along V=0 for Q dep, then use this
%     project-correction thing.. Let's see how that works
%       Ah ok - it just predicts the corresponding w_TO value at (say) V=0
%       and then adds the residual term?
        % which is to say; calculates difference in w_TO between V=v and V=0, then
        % adds this to the actual value

% We could try using the pred errs but they seem very small. But the same
% approach could be taken using the new linear fit.

        
% Alternatively one can use the orthogonal distance and plot the (say) max
% PI along this... 
% %
q_vals = polz_q(use_vals);
v_vals = polz_v(use_vals);
fit_mdl_val = lfit_fin(q_vals,v_vals)*fit_scale_fac+l_fit_offset;
fit_mdl_v_zero = lfit_fin(q_vals,0*v_vals)*fit_scale_fac+l_fit_offset;
fit_mdl_q_m1 = lfit_fin(0*q_vals-1,v_vals)*fit_scale_fac+l_fit_offset;

% 
% opt_theta_fit_call = @(x) opt_theta_eval(polz_cont,polz_theta,polz_v,tune_out_vals_scaled,to_unc,x);
% [theta_L,cost_fun,exitflag,output] = fminunc(opt_theta_fit_call,0);
% [~,opt_out] = opt_theta_eval(polz_cont,polz_theta,polz_v,tune_out_vals_scaled,to_unc,theta_L);
% % [~,opt_out] = opt_theta_eval(polz_cont,polz_theta,polz_v,to_val_fit,to_unc,theta_L)
% % cli_header(1,'Optimized value %.3f',theta_L);
% lfit = opt_out.fit;
% gof = opt_out.gof;
% out = opt_out.out;
% q_A = opt_out.predictors(:,1);
% v_A = opt_out.predictors(:,2);
% predictors = opt_out.predictors;
% observations = opt_out.observations;

lin_pred = lfit_fin([-1,0])*fit_scale_fac+l_fit_offset;
% lin_unc = predint(lfit,[-1,0])-lin_pred;


predict_predictor_input=cat(2,q_vals,v_vals,repmat(-fit_vals_full(4),numel(q_vals),1));
% full_mdl_val = predict(fit_mdl_full,predict_predictor_input);
full_resid = to_fit_resid(use_vals);

% Q_rot = Q_fun(polz_cont,pol

unique_q = unique_polz_vq(:,2);
unique_v = unique_polz_vq(:,1);
% to_binned = 

bin_mdl_val = lfit_fin(unique_q,unique_v)*fit_scale_fac+l_fit_offset;
bin_mdl_v_zero = lfit_fin(unique_q,0*unique_v)*fit_scale_fac+l_fit_offset;
bin_mdl_q_zero = lfit_fin(0*unique_q-1,unique_v)*fit_scale_fac+l_fit_offset;

% define the shift operator s.t. shift + obs = interp
shift_vq=fit_mdl_val-[fit_mdl_v_zero,fit_mdl_q_m1];
bin_shift_vq=bin_mdl_val-[bin_mdl_v_zero,bin_mdl_q_zero];
% we expect that...



q_mdlplot = lfit_fin(Q_plot,0*Q_plot)*fit_scale_fac+l_fit_offset;
v_mdlplot = lfit_fin(0*Q_plot-1,V_plot)*fit_scale_fac+l_fit_offset;
% and to get the uncerts,
           
obs_shift_to_v_0=(observations-shift_vq(:,1)); % GHz
obs_shift_to_q_m1=((observations)-shift_vq(:,2)); % GHz
q_shift_residual = obs_shift_to_v_0 - fit_mdl_v_zero ;
v_shift_residual = obs_shift_to_q_m1 - fit_mdl_q_m1;


intlvls = .01:.02:.95;
show_array = zeros(1,17);
c = 1;
for ii = 1:length(intlvls)-2
    int_lvl = intlvls(ii+1);
    storearray = 1;

    d_pi_nonsim = predint(lfit_fin,[q_vals,v_vals],int_lvl,'observation','off');
    d_pi_sim = predint(lfit_fin,[q_vals,v_vals],int_lvl,'observation','on');
    d_ci_nonsim = predint(lfit_fin,[q_vals,v_vals],int_lvl,'functional','off');
    d_ci_sim = predint(lfit_fin,[q_vals,v_vals],int_lvl,'functional','on');

    
    q_pi_nonsim = predint(lfit_fin,[Q_plot;0*V_plot]',int_lvl,'observation','off');
    q_pi_sim = predint(lfit_fin,[Q_plot;0*V_plot]',int_lvl,'observation','on');
    q_ci_nonsim = predint(lfit_fin,[Q_plot;0*V_plot]',int_lvl,'functional','off');
    q_ci_sim = predint(lfit_fin,[Q_plot;0*V_plot]',int_lvl,'functional','on');
    
    v_pi_nonsim = predint(lfit_fin,[0*Q_plot-1;V_plot]',int_lvl,'observation','off');
    v_pi_sim = predint(lfit_fin,[0*Q_plot-1;V_plot]',int_lvl,'observation','on');
    v_ci_nonsim = predint(lfit_fin,[0*Q_plot-1;V_plot]',int_lvl,'functional','off');
    v_ci_sim = predint(lfit_fin,[0*Q_plot-1;V_plot]',int_lvl,'functional','on');
    
    obs_q_pi_nonsim = predint(lfit_fin,[q_vals,0*v_vals],int_lvl,'observation','off');
    obs_v_pi_nonsim = predint(lfit_fin,[0*q_vals-1,v_vals],int_lvl,'observation','off');
    obs_q_pi_sim = predint(lfit_fin,[q_vals,0*v_vals],int_lvl,'observation','on');
    obs_v_pi_sim = predint(lfit_fin,[0*q_vals-1,v_vals],int_lvl,'observation','on');
    obs_q_ci_nonsim = predint(lfit_fin,[q_vals,0*v_vals],int_lvl,'functional','off');
    obs_v_ci_nonsim = predint(lfit_fin,[0*q_vals-1,v_vals],int_lvl,'functional','off');
    obs_q_ci_sim = predint(lfit_fin,[q_vals,0*v_vals],int_lvl,'functional','on');
    obs_v_ci_sim = predint(lfit_fin,[0*q_vals-1,v_vals],int_lvl,'functional','on');
    
    
%     [full_mdl_val,full_mdl_pi_nonsim]=predict(fit_mdl_full,predict_predictor_input,'Alpha',1-int_lvl,'Prediction','observation','Simultaneous','false');
    [~,full_mdl_pi]=predict(fit_mdl_full,predict_predictor_input,'Alpha',1-int_lvl,'Prediction','observation','Simultaneous','true');
%     [~,full_mdl_ci_nonsim]=predict(fit_mdl_full,predict_predictor_input,'Alpha',1-int_lvl,'Prediction','curve','Simultaneous','false');
    [~,full_mdl_ci]=predict(fit_mdl_full,predict_predictor_input,'Alpha',1-int_lvl,'Prediction','curve','Simultaneous','true');

    % check how many runs fall within the intervals



    % if the residual is within the bounds, then the diff with lower will be
    % >0, diff with upper will be <0, hence product is <0.
    % If above/below the upper/lower bound then the diffs are both the same
    % sign
    v_res_in_intervals = [sum(prod((v_shift_residual-(obs_v_pi_sim-fit_mdl_q_m1)),2)<0)/length(observations),...
                        sum(prod((v_shift_residual-(obs_v_pi_nonsim-fit_mdl_q_m1)),2)<0)/length(observations),...                    
                         sum(prod((v_shift_residual-(obs_v_ci_nonsim-fit_mdl_q_m1)),2)<0)/length(observations),...
                         sum(prod((v_shift_residual-(obs_v_ci_sim-fit_mdl_q_m1)),2)<0)/length(observations)];

    q_res_in_intervals = [sum(prod((q_shift_residual-(obs_q_pi_sim-fit_mdl_v_zero)),2)<0)/length(observations),...
                            sum(prod((q_shift_residual-(obs_q_pi_nonsim-fit_mdl_v_zero)),2)<0)/length(observations),...
                            sum(prod((q_shift_residual-(obs_q_ci_nonsim-fit_mdl_v_zero)),2)<0)/length(observations),...
                            sum(prod((q_shift_residual-(obs_q_ci_sim-fit_mdl_v_zero)),2)<0)/length(observations)];

    data_res_in_intervals = [sum(prod((q_shift_residual-(d_pi_sim-fit_mdl_val)),2)<0)/length(observations),...
                        sum(prod((q_shift_residual-(d_pi_nonsim-fit_mdl_val)),2)<0)/length(observations),...
                        sum(prod((q_shift_residual-(d_ci_nonsim-fit_mdl_val)),2)<0)/length(observations),...
                        sum(prod((q_shift_residual-(d_ci_sim-fit_mdl_val)),2)<0)/length(observations)];
                    
    nl_res_in_intervals = [sum(prod((full_resid-(1e6*full_mdl_pi-fit_mdl_val+to_fit_val_offset)),2)<0)/length(observations),...
                            nan,...
                            sum(prod((full_resid-(1e6*full_mdl_ci-fit_mdl_val+to_fit_val_offset)),2)<0)/length(observations),...
                                nan];
                        

                        
    interval_capture_rate = [q_res_in_intervals',v_res_in_intervals',data_res_in_intervals',nl_res_in_intervals'];
    interval_capture_rate = 100*round(interval_capture_rate,4);

    src = find(show_array(:,1) == int_lvl);
    if ~isempty(src)
        c = src;
    else
        c = c + 1;
    end
    show_array(c,:) = [int_lvl,q_res_in_intervals,v_res_in_intervals,data_res_in_intervals,nl_res_in_intervals];
end


% cli_header('Fraction of shots in %.1fp.c. intervals:',100*int_lvl);
interval_types = ["Simultaneous Pred.","Prediction"',"Confidence","Simult. Conf."]';
capture_table = table(interval_types,interval_capture_rate);



% %

interval_colors = inferno(6);
full_colors = viridis(6);
stfig('Interval captures');
clf
hold on


c = 0;
legs = [0,0,0,0];
for colidx = 10:13
    c = c + 1;
    if c == 1
    dataline = plot(show_array(:,1),show_array(:,colidx),'-','Color',interval_colors(6-c,:),...
        'LineWidth',3);
    else
        plot(show_array(:,1),show_array(:,colidx),'-','Color',interval_colors(6-c,:),...
        'LineWidth',3);
    end
end
c = 0;
for colidx = 2:5
    c = c + 1;
    legs(c) = plot(show_array(:,1),show_array(:,colidx),'o','Color',interval_colors(6-c,:),...
        'MarkerSize',10);
end
c = 0;
for colidx = 6:9
    c = c + 1;
    if c==1
        v_zero_plot = plot(show_array(:,1),show_array(:,colidx),'.','Color',interval_colors(6-c,:),...
            'MarkerSize',20);
    else
        plot(show_array(:,1),show_array(:,colidx),'.','Color',interval_colors(6-c,:),...
        'MarkerSize',20);
    end
end

c = 0;
flegs = [0,0];
for colidx = 14:17
    c = c + 1;
    flegs(c) = plot(show_array(:,1),show_array(:,colidx),'d','Color',full_colors(2+c,:),...
        'MarkerSize',8,'LineWidth',2);
end

plot([0,1],[0,1],':','LineWidth',4,'Color',0.6*[1,1,1])

% legs
legend([dataline,v_zero_plot,legs,flegs],["2D linear model","V=0","Sim. Prediction (Q=0)"',"Prediction","Confidence","Simult. Conf.","Nonlin. pred.","Nonlin. conf."],'Location','best')
xlabel('Confidence level')
ylabel('Capture rate')
title('Fraction of scans with $f_{TO}$ within uncert intervals')
set(gca,'FontSize',20)
box on
set(gca,'linewidth', 1.5)
xticks(0:0.2:1)
yticks(0:0.2:1)




 %polz_vq_to_val_unc columns: Q, V, f_TO, sterr, std
stfig('Q V dep new');
clf
subplot(1,2,1)
hold on
allscan = plot(v_vals,scale_factor_freq*(obs_shift_to_q_m1-pred_TO),'k.','MarkerSize',10);
plot(V_plot,scale_factor_freq*(v_mdlplot-pred_TO),'k')

lp1=patch([V_plot, fliplr(V_plot)],...
    [scale_factor_freq*(v_pi_sim(:,2)-pred_TO)', fliplr(scale_factor_freq*(v_pi_sim(:,1)-pred_TO)')],...
    interval_colors(5,:),'EdgeColor','none','FaceAlpha',0.2);
lp2=patch([V_plot, fliplr(V_plot)],...
    [scale_factor_freq*(v_pi_nonsim(:,2)-pred_TO)', fliplr(scale_factor_freq*(v_pi_nonsim(:,1)-pred_TO)')],...
    interval_colors(4,:),'EdgeColor','none','FaceAlpha',0.2);
lp3=patch([V_plot, fliplr(V_plot)],...
    [scale_factor_freq*(v_ci_sim(:,2)-pred_TO)', fliplr(scale_factor_freq*(v_ci_sim(:,1)-pred_TO)')],...
    interval_colors(3,:),'EdgeColor','none','FaceAlpha',0.2);
lp4=patch([V_plot, fliplr(V_plot)],...
    [scale_factor_freq*(v_ci_nonsim(:,2)-pred_TO)', fliplr(scale_factor_freq*(v_ci_nonsim(:,1)-pred_TO)')],...
    interval_colors(2,:),'EdgeColor','none','FaceAlpha',0.2);
bin_leg = errorbar(polz_vq_to_val_unc(:,1),...
    scale_factor_freq*(polz_vq_to_val_unc(:,3)-pred_TO-bin_shift_vq(:,2)),...
    scale_factor_freq*(polz_vq_to_val_unc(:,5)),...
     'o','CapSize',0,'MarkerSize',5,'Color',colors_main(1,:),...
    'MarkerFaceColor',colors_detail(1,:),'LineWidth',2.8);



hold off
box on
set(gcf,'Units','Pixels')
set(gca,'linewidth', 1.5)
set(gca,'TickLength',[0.02,0])
xlim([min(V_plot),max(V_plot)])
xticks([-1,-0.5,0,0.5,1])
xlabel('$\mathcal{V}$ ($4^{\mathrm{th}}$ Stokes parameter)')
ylabel('$f_\mathrm{TO}(\mathcal{Q_A}=-1,\mathcal{V}) - f_\mathrm{TO}(-1,0)$ (GHz)')
legend([allscan,bin_leg,lp1,lp2,lp3,lp4],{'Scans','Binned data','Sim. PI','PI','Sim CI','CI'},'Location','Best')
set(gca,'FontSize',20)

% 
subplot(1,2,2)
hold on
plot(q_vals,scale_factor_freq*(observations-shift_vq(:,1)-pred_TO),'k.','MarkerSize',10);
plot(Q_plot,scale_factor_freq*(q_mdlplot-pred_TO),'k')

patch([Q_plot, fliplr(Q_plot)],...
    [scale_factor_freq*(q_pi_sim(:,2)-pred_TO)', fliplr(scale_factor_freq*(q_pi_sim(:,1)-pred_TO)')],...
    interval_colors(5,:),'EdgeColor','none','FaceAlpha',0.2)
patch([Q_plot, fliplr(Q_plot)],...
    [scale_factor_freq*(q_pi_nonsim(:,2)-pred_TO)', fliplr(scale_factor_freq*(q_pi_nonsim(:,1)-pred_TO)')],...
    interval_colors(4,:),'EdgeColor','none','FaceAlpha',0.2)
patch([Q_plot, fliplr(Q_plot)],...
    [scale_factor_freq*(q_ci_sim(:,2)-pred_TO)', fliplr(scale_factor_freq*(q_ci_sim(:,1)-pred_TO)')],...
    interval_colors(3,:),'EdgeColor','none','FaceAlpha',0.2)
patch([Q_plot, fliplr(Q_plot)],...
    [scale_factor_freq*(q_ci_nonsim(:,2)-pred_TO)', fliplr(scale_factor_freq*(q_ci_nonsim(:,1)-pred_TO)')],...
    interval_colors(2,:),'EdgeColor','none','FaceAlpha',0.2)


errorbar(polz_vq_to_val_unc(:,2),...
    scale_factor_freq*(polz_vq_to_val_unc(:,3)-pred_TO-bin_shift_vq(:,1)),...
    scale_factor_freq*(polz_vq_to_val_unc(:,5)),...
     scale_factor_freq*(polz_vq_to_val_unc(:,5)),...
     'o','CapSize',0,'MarkerSize',5,'Color',colors_main(1,:),...
    'MarkerFaceColor',colors_detail(1,:),'LineWidth',2.8);

hold off
box on
set(gcf,'Units','Pixels')
set(gca,'linewidth', 1.5)
set(gca,'TickLength',[0.02,0])
xlim([min(Q_plot),max(Q_plot)])
xticks([-1,-0.5,0,0.5,1])
ylabel('$f_\mathrm{TO}(\mathcal{Q_A},\mathcal{V}=0) - f_\mathrm{TO}(-1,0)$ (GHz)')
xlabel('$\mathcal{Q_{A}}$, ($2^{\mathrm{nd}}$ Stokes parameter )') %font_size_label
set(gca,'FontSize',20)


 %% 
 
 
 
conf_level = erf(1/sqrt(2));
 
% q_trend = lfit_lin([Q_plot;0*V_plot]');
q_pi_nonsim = predint(lfit_fin,[Q_plot;0*V_plot]',conf_level,'observation','off');
q_ci_sim = predint(lfit_fin,[Q_plot;0*V_plot]',conf_level,'functional','on');

% v_trend = lfit_lin([0*Q_plot-1;V_plot]');
v_pi_nonsim = predint(lfit_fin,[0*Q_plot-1;V_plot]',conf_level,'observation','off');
v_ci_sim = predint(lfit_fin,[0*Q_plot-1;V_plot]',conf_level,'functional','on');

 
stfig('V dep new');
clf
hold on
% allscan = plot(v_vals,scale_factor_freq*(obs_shift_to_q_m1-pred_TO),'k.','MarkerSize',10);
plot(V_plot,scale_factor_freq*(v_mdlplot-pred_TO),'k','LineWidth',2)
% plot(V_plot
lp2=patch([V_plot, fliplr(V_plot)],...
    [scale_factor_freq*(v_pi_nonsim(:,2)-pred_TO)', fliplr(scale_factor_freq*(v_pi_nonsim(:,1)-pred_TO)')],...
    0.3*[1,1,1],'EdgeColor','none','FaceAlpha',0.2);
% lp3=patch([V_plot, fliplr(V_plot)],...
%     [scale_factor_freq*(v_ci_sim(:,2)-pred_TO)', fliplr(scale_factor_freq*(v_ci_sim(:,1)-pred_TO)')],...
%     0.5*[1,1,1],'EdgeColor','none','FaceAlpha',0.5);
bin_leg = errorbar(polz_vq_to_val_unc(:,1),...
    scale_factor_freq*(polz_vq_to_val_unc(:,3)-pred_TO-bin_shift_vq(:,2)),...
    scale_factor_freq*(polz_vq_to_val_unc(:,5)),...
     'o','CapSize',0,'MarkerSize',5,'Color',colors_main(1,:),...
    'MarkerFaceColor',colors_detail(1,:),'LineWidth',2.8);
pred_pt=errorbar(0,0,(to_scalar_minus_half_tensor.unc+to_scalar_minus_half_tensor.unc_boot)/1e9,...
    '.','MarkerFaceColor','k','MarkerSize',5,'LineWidth',2,'MarkerFaceColor',interval_colors(4,:),...
    'Color',interval_colors(4,:));

hold off
box on
set(gcf,'Units','Pixels')
set(gca,'linewidth', 1.5)
set(gca,'TickLength',[0.02,0])
xlim([min(V_plot),max(V_plot)])
xticks([-1,-0.5,0,0.5,1])
xlabel('$\mathcal{V}$ ($4^{\mathrm{th}}$ Stokes parameter)')
ylabel('$f_\mathrm{TO}(\mathcal{Q_A}=-1,\mathcal{V}) - f_\mathrm{TO}(-1,0)$ (GHz)')
% legend([allscan,bin_leg,pred_pt],{'All scans','Binned data','$f_\mathrm{TO}$(-1,0)'},'Location','SE')
legend([bin_leg,pred_pt],{'Binned data','$f_\mathrm{TO}$(-1,0)'},'Location','SE')
set(gca,'FontSize',20)


stfig('Q dep new');
clf
hold on
% plot(polz_q(use_vals),scale_factor_freq*(observations-shift_vq(:,1)-pred_TO),'k.','MarkerSize',10);
plot(Q_plot,scale_factor_freq*(q_mdlplot-pred_TO),'k','LineWidth',2)


patch([Q_plot, fliplr(Q_plot)],...
    [scale_factor_freq*(q_pi_nonsim(:,2)-pred_TO)', fliplr(scale_factor_freq*(q_pi_nonsim(:,1)-pred_TO)')],...
     0.3*[1,1,1],'EdgeColor','none','FaceAlpha',0.2)
%  patch([Q_plot, fliplr(Q_plot)],...
%     [scale_factor_freq*(q_ci_sim(:,2)-pred_TO)', fliplr(scale_factor_freq*(q_ci_sim(:,1)-pred_TO)')],...
%      0.5*[1,1,1],'EdgeColor','none','FaceAlpha',0.5)

errorbar(polz_vq_to_val_unc(:,2),...
    scale_factor_freq*(polz_vq_to_val_unc(:,3)-pred_TO-bin_shift_vq(:,1)),...
    scale_factor_freq*(polz_vq_to_val_unc(:,5)),...
     scale_factor_freq*(polz_vq_to_val_unc(:,5)),...
     'o','CapSize',0,'MarkerSize',5,'Color',colors_main(1,:),...
    'MarkerFaceColor',colors_detail(1,:),'LineWidth',2.8);

errorbar(-1,0,(to_scalar_minus_half_tensor.unc+to_scalar_minus_half_tensor.unc_boot)/1e9,...
    '.','MarkerFaceColor','k','MarkerSize',5,'LineWidth',2,'MarkerFaceColor',interval_colors(4,:),...
    'Color',interval_colors(4,:))

hold off
box on
set(gcf,'Units','Pixels')
set(gca,'linewidth', 1.5)
set(gca,'TickLength',[0.02,0])
xlim([min(Q_plot),max(Q_plot)])
xticks([-1,-0.5,0,0.5,1])
ylabel('$f_\mathrm{TO}(\mathcal{Q_A},\mathcal{V}=0) - f_\mathrm{TO}(-1,0)$ (GHz)')
xlabel('$\mathcal{Q_{A}}$, ($2^{\mathrm{nd}}$ Stokes parameter )') %font_size_label
set(gca,'FontSize',20)

%% Checking theory values
bT_coef = fit_mdl_full.Coefficients.Estimate(3)*(sin(fit_mdl_full.Coefficients.Estimate(5)))^2;
bV_coef = fit_mdl_full.Coefficients.Estimate(2)*(cos(fit_mdl_full.Coefficients.Estimate(5)));
angle_vals = linspace(20,35,200);
trig_coefs = [sin(deg2rad(angle_vals')).^2,cos(deg2rad(angle_vals'))];
beta_terms = [bT_coef,bV_coef]./trig_coefs;

angle_mean = 27;
angle_sd = 2;
prob_vals = (1/sqrt(2*pi*angle_sd^2)) * exp(-0.5*((angle_vals-angle_mean)/angle_sd).^2);

theta_test = mean(angle_vals);
bT_test = bT_coef/sin(deg2rad(theta_test))^2;
bV_test = bV_coef/cos(deg2rad(theta_test));
bT_theory= 3400;
bV_theory = nan;

min(beta_terms,[],1)
max(beta_terms,[],1)
cli_header('beta^T est %.1f(%.1f,%.1f)',bT_test,min(beta_terms(:,1))-bT_test,max(beta_terms(:,1))-bT_test);
cli_header('beta^V est %.1f(%.1f,%.1f)',bV_test,min(beta_terms(:,2))-bV_test,max(beta_terms(:,2))-bV_test);
cli_header('using angle range %.f-%.f deg',min(angle_vals),max(angle_vals));

stfig('beta');
clf
subplot(2,1,1)
hold on
plot(angle_vals,beta_terms)
plot(theta_test,bT_theory,'kd')
plot(theta_test,bT_test,'kx')
plot(theta_test,bV_test,'ko')

% subplot(2,1,2)
% hold on
% plot(angle_vals,prob_vals.*beta_terms')

function [costfun,opt_out] = opt_theta_eval(polz_cont,polz_theta,v_vals,observations,to_unc,theta_L)

    polz_q=Q_fun(polz_cont,polz_theta,theta_L);
    
     w_fit=1./(to_unc.^2);
     if sum(w_fit)==0
         w_fit = 0*w_fit+1;
     end
     w_fit = w_fit./sum(w_fit);
    
    use_vals = ~isnan(v_vals);
    polz_q = polz_q(use_vals);
    observations = observations(use_vals);
    v_vals = v_vals(use_vals);
    to_unc = to_unc(use_vals);
    w_fit = w_fit(use_vals);
    predictors = [polz_q,v_vals];
    
         
    [lfit,gof,out] = fit(predictors,observations,'poly11','Weights',w_fit);
    
    costfun = gof.rmse;
    opt_out.fit = lfit;
    opt_out.gof = gof;
    opt_out.out = out;
    opt_out.predictors = predictors;
    opt_out.observations = observations;
end

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

predict_predictor_input=cat(2,q_in*0,v_in,repmat(-mdl.Coefficients.Estimate(4),numel(v_in),1));
fit_mdl_q_zero=predict(mdl,predict_predictor_input);

shift_vq=[fit_mdl_v_zero-fit_mdl_val,fit_mdl_q_zero-fit_mdl_val];

end
% 
% function shift_vq=v_lin_shift(v_in,q_in,mdl)
% 
% q_in=col_vec(q_in);
% v_in=col_vec(v_in);
% if numel(v_in)~=numel(q_in)
%     error('v,q must be same size')
% end
% 
% 
% fit_mdl_val=mdl(q_in,v_in);
% 
% fit_mdl_v_zero = mdl(q_in,0*v_in);
% fit_mdl_q_zero = mdl(0*q_in,v_in);
% 
% shift_vq=[fit_mdl_v_zero-fit_mdl_val,fit_mdl_q_zero-fit_mdl_val];
% 
% end

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


use_vals = ~isnan(polz_v);
w_fit_postmask=1./(to_unc(use_vals).^2);
if ~use_fit_weights
    w_fit_postmask=w_fit_postmask*0+1;
end
w_fit_postmask = w_fit_postmask./sum(w_fit_postmask);

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
% update to get full output... 
[beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(predictor,tune_out_vals_scaled,@full_tune_out_polz_model,beta0,...
   'Weights',w_fit,'Options',opts); % identical output from the model


detailed_out.fit=fit_mdl_full;
detailed_out.fit_output.params = beta;
detailed_out.fit_output.residuals = R; 
detailed_out.fit_output.Jacobian= J;
detailed_out.fit_output.MSE = MSE;
detailed_out.fit_output.ErrorModelInfo = ErrorModelInfo;
detailed_out.fit_output.Cov = CovB;


% Trying with rescaled vars?
% Conclusion: Not it. Still ill-conditioned and same predictions.
% also compared to modified model trying to reduce overparametrizasion
% four-parameter fit converges and is not, notably, overparametrized
% i.e. no ill-conditioning warning
% and gives same fit i.e. same residuals i.e same plane

nanmask_keep = ~any(isnan(predictor),2);
detailed_out.nanmask = nanmask_keep;

fit_vals_full = fit_mdl_full.Coefficients.Estimate;
fit_uncs_full = fit_mdl_full.Coefficients.SE;
%
% Trying to eliminate a fit parameter...
% using the following conversions
% a = fTOScalar + \[Beta]T/2 - 3/4 \[Beta]T Sin[\[Theta]k]^2;
% b = 3/4  \[Beta]T Sin[\[Theta]k]^2;
% % c = 1/2  \[Beta]V Cos[\[Theta]k];
% a_reduced_coef = fit_vals_full(1) + fit_vals_full(3)/2 - 0.75*fit_vals_full(3);
% b_reduced_coef = 0.75*fit_vals_full(3)*sin(fit_vals_full(5)).^2;
% c_reduced_coef = 0.5*fit_vals_full(2)*cos(fit_vals_full(5));
% beta_reduced = [a_reduced_coef,b_reduced_coef,c_reduced_coef]'; 
% 
% to_lin_offset = mean(to_val_fit);
% options = optimoptions('fminunc','Display','iter-detailed');
% opt_theta_fit_call = @(x) opt_theta_eval(polz_cont,polz_theta,polz_v,tune_out_vals_scaled,to_unc,x);
% [theta_L,cost_fun,exitflag,output] = fminunc(opt_theta_fit_call,0);
% [~,opt_out] = opt_theta_eval(polz_cont,polz_theta,polz_v,tune_out_vals_scaled,to_unc,theta_L);
% % [~,opt_out] = opt_theta_eval(polz_cont,polz_theta,polz_v,to_val_fit,to_unc,theta_L)
% % cli_header(1,'Optimized value %.3f',theta_L);
% 
% 
% lfit = opt_out.fit;
% gof = opt_out.gof;
% out = opt_out.out;
% % q_A_nl = Q_fun(polz_cont,polz_theta, fit_mdl_full.Coefficients.Estimate(4));
% q_A = opt_out.predictors(:,1);
% v_A = opt_out.predictors(:,2);
% predictors = opt_out.predictors;
% observations = opt_out.observations;
% 
% lin_pred = lfit([-1,0])*1e6 + fit_offset;
% lin_unc = predint(lfit,[-1,0])-lin_pred;




if fit_vals_full(3)<0
    error('fit reduced tensor term is less than zero')
end


% % compare the predictions
% [tsmht_details.scaled_val_predict,tsmht_details.scaled_unc_predict]=...
%     predict(fit_mdl_scaled,[-1,0,-fit_vals_scaled(4)],'Alpha',1-erf(1/sqrt(2)),'Prediction','Curve');
[raw_predict,tsmht_details.unc_predict]=...
    predict(fit_mdl_full,[-1,0,-fit_vals_full(4)],'Alpha',1-erf(1/sqrt(2)),'Prediction','Curve');
    
% There is a problem here - in that the fit value for 'phase' is used as
% the input predictor; however this new model couples this to theta_k as
% well... 

tsmht_details.val_predict=raw_predict*1e6+fit_offset;
tsmht_details.unc_predict=1/2*range(tsmht_details.unc_predict)*1e6;

% tsmht_details.lin_val_predict=lin_pred;
% tsmht_details.lin_unc_predict=lin_unc;

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