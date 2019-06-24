function out_struct=chrip_sine_harmonics_model(tdat,xdat,terms,freq_lims,verbose)
xdat=xdat(:);
tdat=tdat(:);
%given a periodic waveform with lots of harmonics return 
% - phase,amplitude,phase of those harmonics
% - a model
% this code calulates [harmoncis,freq chirp,amp_chirp] terms with a chirped amplitude model for each harmonic

%TODO
% -documentation



if ~isequal(size(freq_lims),[1,2])
    error('freq lims not the right size')
end

options.components_min_amp=1e-4;
options.freq_limits=freq_lims;
 options.components_diff_freq=5;
[fft_pks_sub,dom_freq_det]=dominant_freq_components(tdat,xdat,options);


%TODO: option to allow f_fund*1/2*N
% now see what componets are harmonics of the largest component
fund_freq=fft_pks_sub.freq(1);
multiple_of_fund=fft_pks_sub.freq(1:end)/fund_freq;
is_harmonic=multiple_of_fund>1.6; %remove low freq subharmonics
is_harmonic(1)=true; %allow the fundamental
rounded_harmonic=round(multiple_of_fund);
is_harmonic=is_harmonic & fund_freq*abs(rounded_harmonic-multiple_of_fund)<0.2; %difference in hz
fft_pks_sub.harm_rounded=rounded_harmonic;
%mask  and the fft_pks strucutre


harmonic_terms=terms(1);
%mask out what we will fit
is_harmonic=is_harmonic & (1:numel(is_harmonic))'<=harmonic_terms;
fft_pks_to_fit=struct_mask(fft_pks_sub,is_harmonic);
fft_pks_not_fit=struct_mask(fft_pks_sub,~is_harmonic);


%plot the waveform and the fft and the identified peaks
if verbose>2
    stfig('ac waveform fft','add_stack',1);
    clf
    set(gcf,'color','w')
    subplot(2,1,1)
    plot(tdat,xdat,'b')
    ylabel('AC Mains')
    xlabel('time (s)')
    pause(1e-6)
    subplot(2,1,2)
    semilogy(dom_freq_det.fft_dat(1,:),abs(dom_freq_det.fft_dat(2,:)))
    xlim(freq_lims)
    hold on
    plot(fft_pks_to_fit.freq,fft_pks_to_fit.amp,'gx')
    plot(fft_pks_not_fit.freq,fft_pks_not_fit.amp,'rx')
    hold off
    xlabel('freq (Hz)')
    ylabel('Amplitude (v)')
    legend('components inc. in fit','not included peaks')
end


if sum(diff(sort(fft_pks_to_fit.harm_rounded))==0)>0
    sort(fft_pks_to_fit.harm_rounded)
    error('nonunique harmonics found')
end


%now there are two approaches that could be employed here
% just try to fit with all the harmonics up to order in a single go
% the second would be to build up, fitting to the first and then layering on harmonics from there

% test_param=[0,50,1,0,1,0];
% test_harm=[1,3];
% plot(tdat,harmonic_sine_waves(tdat,test_param,test_harm))


freq_chirp_terms=terms(2);
if numel(terms)<3
    amp_chirp_terms=1; 
else
    amp_chirp_terms=terms(3);
end
% test_param=[0,50,1,0,1,0];
% test_harm=[1,3];
% plot(tdat,harmonic_sine_waves(tdat,test_param,test_harm))
meanx=dom_freq_det.mean_xdat;
stdx=dom_freq_det.std_xdat;



%% do a global fit to a small fraction of the data
% before we call fitnlm we will try to get close to the desired fit parameters using a robust global
% optimizer of a freq-amp chirped sine wave
simple_fit_fun=@(param,time) harmonic_sine_waves_amp_freq_chrip(time,param,fft_pks_to_fit.harm_rounded(1),freq_chirp_terms,amp_chirp_terms);


sub_tlim=[-inf,10/fft_pks_to_fit.freq(1)]; %select the first 10 cycles
idx_sub=fast_sorted_mask(tdat,sub_tlim(1),sub_tlim(2));
xdat_sub=xdat(idx_sub(1):idx_sub(2));
tdat_sub=tdat(idx_sub(1):idx_sub(2));
xstd_sub=std(xdat_sub);

options.components_min_amp=xstd_sub*5e-2;
options.freq_limits=freq_lims;
 options.components_diff_freq=5;
 options.num_components=1;
fft_pks_sub=dominant_freq_components(tdat,xdat,options);

beta0 = [mean(xdat_sub),fft_pks_sub.freq(1),zeros(1,freq_chirp_terms-1)]; %set offset and the terms of the freq taylor series to zero
param_tmp=[fft_pks_sub.phase(1),fft_pks_sub.amp(1)];
%then add all higher order amp chirp terms to be zero
param_tmp=cat(2,param_tmp,zeros(1,amp_chirp_terms-1))';
beta0=[beta0,param_tmp(:)'];

gf_opt=[];
domain_offset_freq_tmp=cat(1,[-1,1]*stdx,... %amp
                [45,55]);         %freq
domain_freq_chirp_tmp=repmat([-1,1],freq_chirp_terms-1,1);
domain_phase_tmp=[-2,2]*pi;
domain_amp_chirp_tmp=cat(1,[0.1,10]*fft_pks_to_fit.amp(1),...
                        repmat([-1,1]*fft_pks_to_fit.amp(1),amp_chirp_terms-1,1));     
gf_opt.domain=cat(1,domain_offset_freq_tmp,...
                    domain_freq_chirp_tmp,...
                    domain_phase_tmp,...
                    domain_amp_chirp_tmp);
gf_opt.start=beta0;
gf_opt.rmse_thresh=stdx*5e-2;
gf_opt.plot=false;
gf_opt.level=2;
gf_out=global_fit(tdat_sub,xdat_sub,simple_fit_fun,gf_opt);

%% global fit with more cycles
% TODO: this process of progressive fitting could be made into its own
% function
sub_tlim=[-inf,100/fft_pks_to_fit.freq(1)]; %select the first 10 cycles
idx_sub=fast_sorted_mask(tdat,sub_tlim(1),sub_tlim(2));
xdat_sub=xdat(idx_sub(1):idx_sub(2));
tdat_sub=tdat(idx_sub(1):idx_sub(2));
gf_opt.start=gf_out.params;
gf_out=global_fit(tdat_sub,xdat_sub,simple_fit_fun,gf_opt);


%% set up for the full fit

fit_fun=@(param,time) harmonic_sine_waves_amp_freq_chrip(time,param,fft_pks_to_fit.harm_rounded,freq_chirp_terms,amp_chirp_terms);
beta0 = [gf_out.params(1),gf_out.params(2:1+freq_chirp_terms)]; %set offset and the terms of the freq taylor series to zero
param_tmp=[fft_pks_to_fit.phase,fft_pks_to_fit.amp];
%then add all higher order amp chirp terms to be zero
param_tmp=cat(2,param_tmp,zeros(numel(fft_pks_to_fit.phase),amp_chirp_terms-1));
param_tmp(1,:)=gf_out.params(end-amp_chirp_terms:end);
param_tmp=param_tmp';
beta0=[beta0,param_tmp(:)'];


%generate the parameter names, 
coef_names={'offset '};
freq_names=arrayfun(@(harm) sprintf('f_d%u',harm),0:(freq_chirp_terms-1),'UniformOutput',false);
coef_names=[coef_names,freq_names];
phase_names=arrayfun(@(harm) sprintf('phase_h%u',harm),fft_pks_to_fit.harm_rounded,'UniformOutput',false);

make_amp_term_labels=@(harm) arrayfun(@(amp_deriv) sprintf('amp_h%u_d%u',harm,amp_deriv),...
                                0:(amp_chirp_terms-1),'UniformOutput',false);

amp_names=arrayfun((make_amp_term_labels),fft_pks_to_fit.harm_rounded,'UniformOutput',false)';


nametmp=cat(2,phase_names,cat(1,amp_names{:}))';
coef_names=cat(2,coef_names,nametmp(:)');

opts = statset('nlinfit');
%opts.MaxIter=0; %use for debuging the inital guess
fit_mdl = fitnlm(tdat,xdat,fit_fun,beta0,'Options',opts,'CoefficientNames',coef_names);
    
%% undo the folding of parameters into the fit function and output as a struct

fit_terms_est=fit_mdl.Coefficients.Estimate;
fit_terms_se=fit_mdl.Coefficients.SE;
fit_terms_out=[];
fit_terms_out.offset.Estimate=fit_terms_est(1);
fit_terms_out.offset.SE=fit_terms_se(1);
fit_terms_out.freq_taylor.Estimate=fit_terms_est(2:freq_chirp_terms+1);
fit_terms_out.freq_taylor.SE=fit_terms_se(2:freq_chirp_terms+1);

phase_amps_est=fit_terms_est(freq_chirp_terms+2:end);
phase_amps_est=reshape(phase_amps_est,amp_chirp_terms+1,[])';

phase_amps_se=fit_terms_se(freq_chirp_terms+2:end);
phase_amps_se=reshape(phase_amps_se,amp_chirp_terms+1,[])';

fit_terms_out.component_phase.Estimate=phase_amps_est(:,1);
fit_terms_out.component_phase.SE=phase_amps_se(:,1);

fit_terms_out.component_amp_taylor.Estimate=phase_amps_est(:,2:end);
fit_terms_out.component_amp_taylor.SE=phase_amps_se(:,2:end);

fit_terms_out.component_harm=fft_pks_to_fit.harm_rounded;

%build the fit function with the best fit values

harm_freq_amp_chirp=cat(2,fit_terms_out.component_harm,...
                        fit_terms_out.component_phase.Estimate,...
                        fit_terms_out.component_amp_taylor.Estimate);
models_out.all_terms=@(x) sum_sine_waves_amp_freq_chirped(x,harm_freq_amp_chirp,...
                            fit_terms_out.offset.Estimate,...
                            fit_terms_out.freq_taylor.Estimate);                      
models_out.fundamental=@(x) sum_sine_waves_amp_freq_chirped(x,harm_freq_amp_chirp(1,:),...
                            fit_terms_out.offset.Estimate,...
                            fit_terms_out.freq_taylor.Estimate);
models_out.var_terms=@(x,n) sum_sine_waves_amp_freq_chirped(x,harm_freq_amp_chirp(1:n,:),...
                            fit_terms_out.offset.Estimate,...
                            fit_terms_out.freq_taylor.Estimate);                           
%models_out.fit_mdl=fit_mdl;    %this is pretty large ~3.5mb for 4s of ac waveform so we will not include                     
         
fit_perf=[];
fit_perf.MSE=fit_mdl.MSE;
fit_perf.RMSE=fit_mdl.RMSE;
fit_perf.NumObservations=fit_mdl.NumObservations;
fit_perf.Rsquared=fit_mdl.Rsquared;


out_struct.models=models_out;
out_struct.terms=fit_terms_out;
out_struct.fit_perf=fit_perf;


if verbose>2
    %%plot the result
    %sfigure(3);
    stfig('ac waveform fit','add_stack',1);
    clf
    subplot(2,1,1)
    plot(tdat,xdat,'k')
    hold on
    [y_fit_val,y_ci_fit]=predict(fit_mdl,tdat,'Prediction' ,'observation');
    plot(tdat,y_fit_val,'r')
    plot(tdat,y_ci_fit,'b')
    plot(tdat,models_out.fundamental(tdat),'m')
    plot(tdat,models_out.var_terms(tdat,3),'c')
    plot(tdat,models_out.all_terms(tdat),'g')
    first_guess=fit_fun(beta0,tdat);
    plot(tdat,first_guess,'g')
    legend('data','fit','fit ci','fundemental model','3 freq model','all freq model','first guess')
    hold off
    xlim([0,0.1])
    ylabel('voltage');
    xlabel('time');

    subplot(2,1,2)
    plot(tdat,xdat-y_fit_val,'k')
    ylabel('residuals');
    xlabel('time');
end


end



function out=harmonic_sine_waves_amp_freq_chrip(x,param,harmonic,num_freq_terms,num_amp_terms)
%param = offset,fund_freq_d0,fund_freq_d1,...fund_freq_d(num_freq_terms-1),...
%phase1,amp1_d0,amp1_d1,...,amp1_d(num_amp_terms-1),...
%phase2,amp2_d0,amp2_d1,...,amp2_d(num_amp_terms-1),...
%phaseN,ampN_d0,...,ampN_d(num_amp_terms-1)
% where d0,d1,dj indicates jth derivative
if num_freq_terms==0
    error('cant have no freq terms')
end
num_params=numel(param);
if mod(num_params-1-num_freq_terms,1+num_amp_terms)~=0
    error('params length must be 1+num_freq_terms+N*(1+num_amp_terms)')
end
if fix((num_params-1-num_freq_terms)/(1+num_amp_terms))~=numel(harmonic)
    error('number of elements in harmonic vector wrong')
end

%strip out the parts of the parameters
offset=param(1);
param=param(2:end); %strip this off
freq_terms=param(1:num_freq_terms);
param=param(num_freq_terms+1:end); %the phase,amp tay series

phase_amp=reshape(param,1+num_amp_terms,[])';
harm_freq_amp_chirp=[harmonic(:),phase_amp(:,1),phase_amp(:,2:end)];

out=sum_sine_waves_amp_freq_chirped(x,harm_freq_amp_chirp,offset,freq_terms);

end

