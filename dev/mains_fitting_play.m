% mains waveform fit play
%trying to see what can be done in terms of fitting a function to the recorded mains waveform

load('ai_log_state_before_is_sm.mat')




%%

ac_mains_dat=ai_dat.Data(5,:);
ac_mains_data_filt=gaussfilt(ai_dat.time,ac_mains_dat,5e-4);

ac_mains_mean=mean(ac_mains_dat);
ac_mains_std=std(ac_mains_dat);
%% build a function to automate the process below
fft_plot_freq_lims=[5,2000];
make_harmonics_model(ai_dat.time,ac_mains_dat,4,fft_plot_freq_lims,5)


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