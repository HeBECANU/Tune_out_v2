%tune_out_polz_dep_unc_weight
% this needs to be run after tune_out_polz_dep

% i want to look what happens if we weight by the uncertianty in a data
% run using the polarization state uncert


%% the run data
 
to_val_for_polz=data.main.lin.to.val;
to_unc=data.main.lin.to.unc;
wlin=1./(to_unc.^2);
wlin_fit = wlin./sum(wlin);
to_fit_val_offset=wmean(to_val_for_polz,wlin_fit); %should be very insensitive to this


%% find the polarization for each run
%polarisation model/data options
pol_opts.location = 'post';%post, pre_cen, pre_left, pre_right
pol_opts.predict_method = 'only_data';%'full_fit_pref_data','full_fit_pref_data','full_fit_only','only_data'; %obs (obsovation) fit (pertial fit) full_fit (fit with all parameters free)
                                        %'interp_only','interp_pref_data','interp_pref_interp'
                                        %'gauss_only','gauss_pref_data','gauss_pref_interp'
pol_opts.smoothing=3; %deg
pol_opts.wrap_hwp=0;

pol_opts.hwp=data.main.wp.hwp;
pol_opts.qwp=data.main.wp.qwp;
% 
% pol_opts.hwp=data.main.wp.hwp;
% pol_opts.qwp=data.main.wp.qwp;

pol_opts.add_noise=0; % used as a numerical way of propagating uncert
pol_model=pol_data_query(pol_opts);
%pol_model_basic=pol_data_query_basic(pol_opts);


%% use a basic numerical test of the polarization state CI
% this checks the propagation of polarization measurmeent by adding the esitmated error to the
% polarizer measuremetns

do_polz_prop_numericaly=1

if do_polz_prop_numericaly
    num_repeats=30;
    noised_pol_models=cell(num_repeats,1);
    pol_opts.add_noise=1; % used as a numerical way of propagating uncert
    for ii=1:num_repeats
        
        noised_pol_models{ii}=pol_data_query(pol_opts);
    end
    
    
    noised_v=cellfun(@(x) x.v.val,noised_pol_models,'UniformOutput',0);
    noised_v=cat(2,noised_v{:});
    noised_v_cis=prctile(noised_v',[1-erf(1/sqrt(2)),erf(1/sqrt(2))]*100);
    noised_v_cis=noised_v_cis';
    noised_v_cis=noised_v_cis-repmat(mean(noised_v,2),[1,2])
    
    noised_cont=cellfun(@(x) x.cont.val,noised_pol_models,'UniformOutput',0);
    noised_cont=cat(2,noised_cont{:});
    noised_cont_cis=prctile(noised_cont',[1-erf(1/sqrt(2)),erf(1/sqrt(2))]*100);
    noised_cont_cis=noised_cont_cis';
    noised_cont_cis=noised_cont_cis-repmat(mean(noised_cont,2),[1,2]);
    
    noised_v_std=std(noised_v,[],2);
    stfig('ci eval')
    subplot(2,2,1)
    plot(pol_model.v.val,pol_model.v.ci(:,1)./noised_v_cis(:,1),'x')
    xlabel('v value')
    ylabel('ci basic prop/ num exp  (-ve ci)')
    subplot(2,2,2)
    plot(pol_model.v.val,pol_model.v.ci(:,2)./noised_v_cis(:,2),'x')
    xlabel('v value')
    ylabel('ci basic prop/ num exp  (+ve ci)')
    
    subplot(2,2,3)
    plot(pol_model.v.val,pol_model.cont.ci(:,1)./noised_cont_cis(:,1),'x')
    xlabel('cont value')
    ylabel('ci basic prop/ num exp  (-ve ci)')
    subplot(2,2,4)
    plot(pol_model.v.val,pol_model.cont.ci(:,2)./noised_cont_cis(:,2),'x')
    xlabel('cont value')
    ylabel('ci basic prop/ num exp  (+ve ci)')
    
    % set what we have found as the polarization state error
    pol_model.cont.ci=noised_cont_cis;
    pol_model.v.ci=noised_v_cis;
else
    % to save time we can correct the rough propagation by dividing by the rough ratio
    % dont use this for the final plots
    warning('for prototype only')
    pol_model.cont.ci=pol_model.cont.ci/2;
    pol_model.v.ci=pol_model.v.ci/2;
end


%% combine polz and stat uncert
% use the previous fit to propagate the polarization uncert into freq uncert

for ii=1:3
% Find the polarization q paramter in the atomic frame
%we will compute the second Stokes parameter for the data
pol_model.q_at=[];
[pol_model.q_at.val,pol_model.q_at.ci]=polz_Q_fun_unc(pol_model.cont.val,pol_model.theta.val,fit_param_vals_full(4),...
                                pol_model.cont.ci,pol_model.theta.ci,fit_param_vals_unc(4));


deriv_to_with_v=derivest(@(x) predict(fit_mdl_full,[0,x,-fit_param_vals_full(4)]),1,'Vectorized','no');
deriv_to_with_v=deriv_to_with_v*1e6;

deriv_to_with_q=derivest(@(x) predict(fit_mdl_full,[x,0,-fit_param_vals_full(4)]),1,'Vectorized','no');
deriv_to_with_q=deriv_to_with_q*1e6;
deriv_to_with_q=abs(deriv_to_with_q);

% compine the measurment error with the propagated error in the v state
combined_to_ci_vqdep=(pol_model.v.ci*deriv_to_with_v.*[-1,1]).^2+(fliplr(pol_model.q_at.ci.*[-1,1])*deriv_to_with_q).^2+([-to_unc,to_unc]).^2;
combined_to_ci_vqdep=[-sqrt(combined_to_ci_vqdep(:,1)),sqrt(combined_to_ci_vqdep(:,2))];

% convert this ci into just sysmteric std
combined_to_unc_vqdep=rms(combined_to_ci_vqdep,2);


% do a fit

to_val_fit = to_val_for_polz;
wlin_fit = wlin./sum(wlin);
to_fit_val_offset=wmean(to_val_fit,wlin_fit); %should be very insensitive to this

fit_order=1;
use_weights=1;
% set up for the bootstrap wrapper
full_data_set=[pol_model.cont.val,pol_model.v.val,pol_model.theta.val,to_val_for_polz,combined_to_unc_vqdep];
full_data_set=full_data_set(~any(isnan(full_data_set),2),:)
full_data_set=num2cell(full_data_set,2);
% check that function gives the same answer as the above procedure
[wrap_fun_ans,detailed_out]=two_stage_two_dim_to_fit(full_data_set,to_fit_val_offset,use_weights,fit_order);
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

pause(0)
end

%%

% Do the bootstrap

% %detailed bootstrap
boot=bootstrap_se(@two_stage_two_dim_to_fit,full_data_set,...
    'opp_arguments',{to_fit_val_offset,use_weights,fit_order},...
    'plots',true,...
    'replace',true,...
    'plot_fig_name','TO boot MHz',...
    'samp_frac_lims',[0.10,1],...%[0.005,0.9]
    'num_samp_frac',5,... %20
    'num_samp_rep',10,... %1e2
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
