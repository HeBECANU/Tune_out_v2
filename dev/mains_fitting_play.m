% mains waveform fit play
%trying to see what can be done in terms of fitting a function to the recorded mains waveform

load('ai_log_state_before_is_sm.mat')


%Goals
% - understand how mains harmonics works
% - recover the 



%%

ac_mains_dat=ai_dat.Data(5,:);
ac_mains_data_filt=gaussfilt(ai_dat.time,ac_mains_dat,5e-4);

ac_mains_mean=mean(ac_mains_dat);
ac_mains_std=std(ac_mains_dat);
%% build a function to automate the process below
fft_plot_freq_lims=[5,5000];
make_harmonics_model(ai_dat.time,ac_mains_dat,[1,3],fft_plot_freq_lims,5)

%% try with harmonics,freq chirp and amp chirp
fft_plot_freq_lims=[5,5000];
complicated_harmonics_model(ai_dat.time,ac_mains_dat,[20,2,2],fft_plot_freq_lims,5)
%given 20 harmonics,3 freq terms how does the numer of amp terms change the RMSE
% 1 0.00317
% 2 0.00306
% 3 0.00306

% given 20 harmonics 2 amp terms how do the freq terms change things
% 3 0.00306
% 2 0.0031
% 1 0.00906


%% can we improve things by adding a time dependent amplitude
% how do the harmonic amplidudes relate to the fundemental amplitude

%%

fft_mains=fft_tx(ai_dat.time,ac_mains_dat-ac_mains_mean,'window','blackman','padding',100);

%fig_handle=figure('Name','Ac input','NumberTitle','off');
sfigure(2);
clf
set(gcf,'color','w')
subplot(2,1,1)
plot(ai_dat.time,ac_mains_dat,'b')
hold on
plot(ai_dat.time,ac_mains_data_filt,'k')
ylabel('AC Mains')
xlabel('time (s)')
pause(1e-6)
subplot(2,1,2)
semilogy(fft_mains(1,:),abs(fft_mains(2,:)))
xlim(fft_plot_freq_lims)

ai_log_single_out.ac_model=nan;

% find peaks
[pks_amp,pks_idx] = findpeaks(abs(fft_mains(2,:)),'MinPeakHeight',ac_mains_std/1000,'MinPeakDistance',round(2/diff(fft_mains(1,1:2))));
pks_freq=fft_mains(1,pks_idx);
[~,sort_order]=sort(pks_amp,'descend');
pks_amp=pks_amp(sort_order);
pks_freq=pks_freq(sort_order);
pks_phase=angle(fft_mains(2,pks_idx))+pi/2;
%get the index



hold on
plot(pks_freq,pks_amp,'rx')
hold off


%% instfreq

[ifq,tifq]=instfreq(ac_mains_data_filt, ai_dat.sample_rate);
plot(tifq,ifq)


%% try to fit a sine wave based on the FFT

xdat=ai_dat.time;
ydat=ac_mains_dat;
simple_sine_fun = @(x,amp,freq,phase,offset) amp.*sin(x.*2*pi*freq+phase)+offset;
fit_fun=@(param,x) simple_sine_fun(x,param(1),param(2),param(3),param(4));
opts = statset('nlinfit');
%opts.MaxIter=0;
%opts.RobustWgtFun = 'welsch' ; %a bit of robust fitting
%opts.Tune = 1;
beta0 = [pks_amp(1),pks_freq(1),pks_phase(1),ac_mains_mean]; %intial guesses
fit_mdl = fitnlm(xdat,ydat,fit_fun,beta0,'Options',opts,'CoefficientNames',{'amp','freq','phase','offset'});
x_samp_fit=ai_dat.time;
[y_fit_val,y_ci_fit]=predict(fit_mdl,x_samp_fit');
sfigure(2);
subplot(2,1,1)
hold on
plot(x_samp_fit,y_fit_val,'r')
plot(x_samp_fit,y_ci_fit,'g')
hold off

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



%% try to fit with a few harmonics


xdat=ai_dat.time;
ydat=ac_mains_dat;
simple_sine_fun = @(x,amp,freq,phase,offset,amp3,phase3)...
    amp.*sin(x.*2*pi*freq+phase)+offset+amp3*sin(x.*3*pi*freq*2+phase3);
fit_fun=@(param,x) simple_sine_fun(x,param(1),param(2),param(3),param(4),param(5),param(6));
opts = statset('nlinfit');
%opts.MaxIter=0;
%opts.RobustWgtFun = 'welsch' ; %a bit of robust fitting
%opts.Tune = 1;
beta0 = [pks_amp(1),pks_freq(1),pks_phase(1),ac_mains_mean,pks_amp(2),pks_phase(2)]; %intial guesses
fit_mdl = fitnlm(xdat,ydat,fit_fun,beta0,'Options',opts,'CoefficientNames',{'amp','freq','phase','offset','amp3','phase3'});
x_samp_fit=ai_dat.time;
[y_fit_val,y_ci_fit]=predict(fit_mdl,x_samp_fit');
sfigure(2);
subplot(2,1,1)
hold on
plot(x_samp_fit,y_fit_val,'r')
plot(x_samp_fit,y_ci_fit,'g')
hold off

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


%% try it with a few chirp terms


xdat=ai_dat.time;
ydat=ac_mains_dat;
simple_sine_fun = @(x,amp,freq,phase,offset,chirp1,chirp2)...
    amp.*sin(2*pi.*x.*(freq+x.*chirp1+(x.^2).*chirp2)+phase)+offset;
fit_fun=@(param,x) simple_sine_fun(x,param(1),param(2),param(3),param(4),param(5),param(6));
opts = statset('nlinfit');
%opts.MaxIter=0;
%opts.RobustWgtFun = 'welsch' ; %a bit of robust fitting
%opts.Tune = 1;
beta0 = [pks_amp(1),pks_freq(1),pks_phase(1),ac_mains_mean,0,0]; %intial guesses
fit_mdl = fitnlm(xdat,ydat,fit_fun,beta0,'Options',opts,'CoefficientNames',{'amp','freq','phase','offset','chirp1','chirp2'});
x_samp_fit=ai_dat.time;
[y_fit_val,y_ci_fit]=predict(fit_mdl,x_samp_fit');
sfigure(2);
subplot(2,1,1)
hold on
plot(x_samp_fit,y_fit_val,'r')
plot(x_samp_fit,y_ci_fit,'g')
hold off

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



%%
num_harm_fit=20;

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
[y_fit_val,y_ci_fit]=predict(fit_mdl,tdat,'Prediction' ,'observation');
sfigure(2);
clf
subplot(2,1,1)
plot(tdat,xdat,'k')
hold on
plot(tdat,y_fit_val,'r')
plot(tdat,y_ci_fit,'b')
hold off

subplot(2,1,2)
plot(tdat,xdat-y_fit_val,'k')


print_var=@(idx) sprintf('%s=%.2f±%.2f',...
            fit_mdl.Coefficients.Row{idx},...
            fit_mdl.Coefficients.Estimate(idx),...
            fit_mdl.Coefficients.SE(idx) );


