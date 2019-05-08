function [harmonic_data,model]=make_harmonics_model(tdat,xdat,order,freq_lims,verbose)
xdat=xdat(:);
tdat=tdat(:);
%given a periodic waveform with lots of harmonics return 
% - phase,amplitude,phase of those harmonics
% - a model

min_peak_factor=1e-4; %amplitude of freq component relative to rms to be included

if size(freq_lims)~=[1,2] 
    error('freq lims not the right size')
end

mean_xdat=mean(xdat);
std_xdat=std(xdat);

fft_dat=fft_tx(tdat,xdat-mean_xdat,'window','chebyshev','win_param',{300},'padding',100);

%mask the fft data to lie within freq_lims
fft_idx_lims=fast_sorted_mask(fft_dat(1,:),freq_lims(1),freq_lims(2));
fft_dat=fft_dat(:,fft_idx_lims(1):fft_idx_lims(2));


% find peaks that are min_peak_factor*xstd and seperated by at least a few times the time resolution
min_pk_sep=5/diff(tdat([1,end])); %peak sep in hz
min_pk_sep_idx=round(min_pk_sep/diff(fft_dat(1,1:2))); %peak sep in fft bins
[pks_unsorted,pks_idx] = findpeaks(abs(fft_dat(2,:)),'MinPeakHeight',std_xdat*min_peak_factor,'MinPeakDistance',min_pk_sep_idx);
pks_freq=fft_dat(1,pks_idx);
%amplitude sort the peaks from highest to smallest amplitude
[~,sort_order]=sort(pks_unsorted,'descend');
fft_pks=[];
fft_pks.amp=pks_unsorted(sort_order);
fft_pks.freq=pks_freq(sort_order);
fft_pks.phase=angle(fft_dat(2,pks_idx))+pi/2;


%plot the waveform and the fft and the identified peaks
if verbose>2
    sfigure(2);
    clf
    set(gcf,'color','w')
    subplot(2,1,1)
    plot(tdat,xdat,'b')
    ylabel('AC Mains')
    xlabel('time (s)')
    pause(1e-6)
    subplot(2,1,2)
    semilogy(fft_dat(1,:),abs(fft_dat(2,:)))
    xlim(freq_lims)
    hold on
    plot(fft_pks.freq,fft_pks.amp,'rx')
    hold off
    xlabel('freq (Hz)')
    ylabel('Amplitude (v)')
end

%TODO: option to allow f_fund*1/2*N
% now see what componets are harmonics of the largest component
multiple_of_fund=fft_pks.freq(1:end)/fft_pks.freq(1);
is_harmonic=multiple_of_fund>1.6; %remove low freq subharmonics
is_harmonic(1)=true; %allow the fundemental
rounded_harmonic=round(multiple_of_fund);
is_harmonic=is_harmonic & abs(rounded_harmonic-multiple_of_fund)<0.005;
fft_pks.harm_rounded=rounded_harmonic;
%mask  and the fft_pks strucutre
fft_pks=struct_mask(fft_pks,is_harmonic);

%now there are two approaches that could be employed here
% just try to fit with all the harmonics up to order in a single go
% the second would be to build up, fitting to the first and then layering on harmonics from there

num_harm_fit=3;

test_param=[0,50,1,0,1,0];
test_harm=[1,3];
plot(tdat,harmonic_sine_waves(tdat,test_param,test_harm))

opts = statset('nlinfit');
%opts.MaxIter=0;
%opts.RobustWgtFun = 'welsch' ; %a bit of robust fitting
%opts.Tune = 1;

fit_fun=@(param,time) harmonic_sine_waves(time,param,fft_pks.harm_rounded(1:num_harm_fit));
beta0 = [0,fft_pks.freq(1)]; %intial guesses
param_tmp=[fft_pks.amp(1:num_harm_fit);fft_pks.phase(1:num_harm_fit)];
beta0=[beta0,param_tmp(:)'];

coef_names={'offset ','freq'};
amp_names=arrayfun(@(harm) sprintf('amp%u',harm),rounded_harmonic(1:num_harm_fit),'UniformOutput',false);
phase_names=arrayfun(@(harm) sprintf('phase%u',harm),rounded_harmonic(1:num_harm_fit),'UniformOutput',false);
nametmp=[amp_names;phase_names];
coef_names=[coef_names,nametmp(:)'];

fit_mdl = fitnlm(tdat,xdat,fit_fun,beta0,'Options',opts,'CoefficientNames',coef_names);
[y_fit_val,y_ci_fit]=predict(fit_mdl,tdat');
sfigure(2);
subplot(2,1,1)
plot(tdat,xdat,'k')
hold on
plot(x_samp_fit,y_fit_val,'r')
plot(x_samp_fit,y_ci_fit,'b')
hold off

subplot(2,1,2)
plot(tdat,xdat-y_fit_val,'k')


print_var=@(idx) sprintf('%s=%.2f±%.2f',...
            fit_mdl.Coefficients.Row{idx},...
            fit_mdl.Coefficients.Estimate(idx),...
            fit_mdl.Coefficients.SE(idx) );
        
% add a box with the fit param
dim = [0.6 0.5 0.3 0.3];
str = {print_var(1),...
      print_var(2),...
       print_var(3)};
annotation('textbox',dim,'String',str,'FitBoxToText','on');




end


function out=harmonic_sine_waves(x,param,harmonic)
%param offset,fund_freq,amp1,phase1,amp2,phase2,...,ampN,freqN,phaseN
num_params=numel(param);
if mod(num_params-2,2)~=0
    error('params length must be 2+N*3')
end
if fix((num_params-2)/2)~=numel(harmonic)
    error('number of elements in harmonic vector wrong')
end
amp_phase=reshape(param(3:end),2,[])';
amp_freq_phase=[amp_phase(:,1),harmonic(:)*param(2),amp_phase(:,2)];
out=sum_sine_waves(x,amp_freq_phase,param(1));

end
