function [harmonic_data,model]=make_harmonics_model(tdat,xdat,num_harm_chirp_fit,freq_lims,verbose)
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


%TODO: option to allow f_fund*1/2*N
% now see what componets are harmonics of the largest component
fund_freq=fft_pks.freq(1);
multiple_of_fund=fft_pks.freq(1:end)/fund_freq;
is_harmonic=multiple_of_fund>1.6; %remove low freq subharmonics
is_harmonic(1)=true; %allow the fundemental
rounded_harmonic=round(multiple_of_fund);
is_harmonic=is_harmonic & fund_freq*abs(rounded_harmonic-multiple_of_fund)<0.2; %difference in hz
fft_pks.harm_rounded=rounded_harmonic;
%mask  and the fft_pks strucutre

%mask out what we will fit
is_harmonic=is_harmonic & (1:numel(is_harmonic))<=num_harm_chirp_fit(1);
fft_pks_to_fit=struct_mask(fft_pks,is_harmonic);
fft_pks_not_fit=struct_mask(fft_pks,~is_harmonic);


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

if numel(num_harm_chirp_fit)==1
    opts = statset('nlinfit');
    %opts.MaxIter=0;
    fit_fun=@(param,time) harmonic_sine_waves(time,param,fft_pks_to_fit.harm_rounded);
    beta0 = [0,fft_pks_to_fit.freq(1)]; %intial guesses
    param_tmp=[fft_pks_to_fit.amp;fft_pks_to_fit.phase];
    beta0=[beta0,param_tmp(:)'];
    coef_names={'offset ','freq'};
    amp_names=arrayfun(@(harm) sprintf('amp%u',harm),fft_pks_to_fit.harm_rounded,'UniformOutput',false);
    phase_names=arrayfun(@(harm) sprintf('phase%u',harm),fft_pks_to_fit.harm_rounded,'UniformOutput',false);
    nametmp=[amp_names;phase_names];
    coef_names=[coef_names,nametmp(:)'];
    fit_mdl = fitnlm(tdat,xdat,fit_fun,beta0,'Options',opts,'CoefficientNames',coef_names)
else
    num_freq_terms=num_harm_chirp_fit(2);
    % test_param=[0,50,1,0,1,0];
    % test_harm=[1,3];
    % plot(tdat,harmonic_sine_waves(tdat,test_param,test_harm))
    opts = statset('nlinfit');
    %opts.MaxIter=0;
    fit_fun=@(param,time) harmonic_sine_waves_chrip(time,param,fft_pks_to_fit.harm_rounded,num_freq_terms);
    beta0 = [0,fft_pks_to_fit.freq(1),zeros(1,num_freq_terms-1)]; %set offset and the terms of the freq taylor series to zero
    param_tmp=[fft_pks_to_fit.amp;fft_pks_to_fit.phase];
    beta0=[beta0,param_tmp(:)'];
    coef_names={'offset '};
    freq_names=arrayfun(@(harm) sprintf('f^(%u)',harm),0:(num_freq_terms-1),'UniformOutput',false);
    coef_names=[coef_names,freq_names];
    amp_names=arrayfun(@(harm) sprintf('amp%u',harm),fft_pks_to_fit.harm_rounded,'UniformOutput',false);
    phase_names=arrayfun(@(harm) sprintf('phase%u',harm),fft_pks_to_fit.harm_rounded,'UniformOutput',false);
    nametmp=[amp_names;phase_names];
    coef_names=[coef_names,nametmp(:)'];
    fit_mdl = fitnlm(tdat,xdat,fit_fun,beta0,'Options',opts,'CoefficientNames',coef_names)
    
end


[y_fit_val,y_ci_fit]=predict(fit_mdl,tdat,'Prediction' ,'observation');
sfigure(3);
clf
subplot(2,1,1)
plot(tdat,xdat,'k')
hold on
plot(tdat,y_fit_val,'r')
plot(tdat,y_ci_fit,'b')
hold off

subplot(2,1,2)
plot(tdat,xdat-y_fit_val,'k')





end


function out=harmonic_sine_waves(x,param,harmonic)
%param offset,fund_freq,amp1,phase1,amp2,phase2,...,ampN,freqN,phaseN
num_params=numel(param);
if mod(num_params-2,2)~=0
    error('params length must be 2+N*2')
end
if fix((num_params-2)/2)~=numel(harmonic)
    error('number of elements in harmonic vector wrong')
end
amp_phase=reshape(param(3:end),2,[])';
amp_freq_phase=[amp_phase(:,1),harmonic(:)*param(2),amp_phase(:,2)];
out=sum_sine_waves(x,amp_freq_phase,param(1));

end

function out=harmonic_sine_waves_chrip(x,param,harmonic,num_freq_terms)
%param offset,fund_freq,amp1,phase1,amp2,phase2,...,ampN,freqN,phaseN
if num_freq_terms==0
    error('cant have no freq terms')
end
num_params=numel(param);
if mod(num_params-1-num_freq_terms,2)~=0
    error('params length must be 1+num_freq_terms+N*2')
end
if fix((num_params-1-num_freq_terms)/2)~=numel(harmonic)
    error('number of elements in harmonic vector wrong')
end

%strip out the parts of the parameters
offset=param(1);
param=param(2:end);
freq_terms=param(1:num_freq_terms);
param=param(num_freq_terms+1:end);

amp_phase=reshape(param,2,[])';
amp_freq_phase=[amp_phase(:,1),harmonic(:),amp_phase(:,2)];
out=sum_sine_waves_chirped(x,amp_freq_phase,offset,freq_terms);

end

