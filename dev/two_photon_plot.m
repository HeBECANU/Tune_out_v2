%make some publication quality plots of the calibration to the 2p transitions

close('all')
%add all subfolders
folder = fileparts(which(mfilename));
folder=strsplit(folder,filesep); %go up a directory
folder=strjoin(folder(1:end-1),filesep);
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));


%%
anal_opt.dir='Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190224_qwp_142\wm_cal';
loz3fun_offset_grad = @(b,x) b(1).*((b(2)/2)^2)./((x-b(3)).^2+(b(2)/2)^2) + b(4)+b(5).*x; %lorentzian

%%
data_dir=dir(fullfile(anal_opt.dir,'2p_cs_log*'));
data_file=data_dir(1).name;
fid=fopen(fullfile(anal_opt.dir,data_file),'r');
raw_file = textscan(fid, '%s','delimiter', '\n');
fclose(fid);
data=jsondecode( raw_file{1}{5});

transition_delta=data.parameters.set_freq-data.parameters.transition_freq;
pmt_current=-1e9*data.parameters.pmt_voltage_mean/data.parameters.current_gain;

opts = statset('MaxIter',1e3);
%opts.RobustWgtFun = 'welsch' ; %a bit of robust fitting
%opts.Tune = 1;
cen_est=wmean(transition_delta,pmt_current);
std_est=std(transition_delta,pmt_current);
amp_est=range(pmt_current);
offset_est=mean(pmt_current)-amp_est/2;
beta0 = [amp_est,std_est,cen_est,offset_est,0]; %intial guesses
fit_mean = fitnlm(transition_delta,pmt_current,loz3fun_offset_grad,beta0,'Options',opts,...
    'CoefficientNames',{'amp','FWHM','center','offset','grad'});
x_samp=linspace(min(transition_delta),max(transition_delta),1e3)';
ci_size_disp=1-erf(1/sqrt(2));%one sd %confidence interval to display
[y_lin,yci_lin]=predict(fit_mean,x_samp,'Prediction','observation','Alpha',ci_size_disp); %'Prediction','observation'
y_mdl=predict(fit_mean,x_samp);

plot(transition_delta,pmt_current-fit_mean.Coefficients.Estimate(4),'.','MarkerEdgeColor',[0.8500, 0.3250, 0.0980])
xlabel('Transition Detuning MHz')
ylabel('Photomultipler current (nA)')
set(gcf,'color','w')

hold on
plot( x_samp,y_mdl-fit_mean.Coefficients.Estimate(4))
hold off

%%
hold on
data_file=data_dir(2).name;
fid=fopen(fullfile(anal_opt.dir,data_file),'r');
raw_file = textscan(fid, '%s','delimiter', '\n');
fclose(fid);
data=jsondecode( raw_file{1}{5});

transition_delta=data.parameters.set_freq-data.parameters.transition_freq;
pmt_current=-1e9*data.parameters.pmt_voltage_mean/data.parameters.current_gain;

opts = statset('MaxIter',1e3);
%opts.RobustWgtFun = 'welsch' ; %a bit of robust fitting
%opts.Tune = 1;
cen_est=wmean(transition_delta,pmt_current);
std_est=std(transition_delta,pmt_current);
amp_est=range(pmt_current);
offset_est=mean(pmt_current)-amp_est/2;
beta0 = [amp_est,std_est,cen_est,offset_est,0]; %intial guesses
fit_mean = fitnlm(transition_delta,pmt_current,loz3fun_offset_grad,beta0,'Options',opts,...
    'CoefficientNames',{'amp','FWHM','center','offset','grad'});
x_samp=linspace(min(transition_delta),max(transition_delta),1e3)';
ci_size_disp=1-erf(1/sqrt(2));%one sd %confidence interval to display
[y_lin,yci_lin]=predict(fit_mean,x_samp,'Prediction','observation','Alpha',ci_size_disp); %'Prediction','observation'
y_mdl=predict(fit_mean,x_samp);

plot(transition_delta,pmt_current-fit_mean.Coefficients.Estimate(4),'k.')
xlabel('Transition Detuning MHz')
ylabel('Photomultipler current (nA)')
set(gcf,'color','w')


hold on
plot( x_samp,y_mdl-fit_mean.Coefficients.Estimate(4))
hold off
