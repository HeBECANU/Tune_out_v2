%try and find the optimal amplitude to sample a sine 

freq_truth =1;
phase_truth =0;
offset_truth=0;
grad_truth  =0;
noise=1;
repeat_fit=10;
amp_list=linspace(0.01,20,1e3);

samples=1e2;
time_max=1.5;
noise_amp=1;
fit_plots=false;

modelfun=@(b,x) b(1).*0.5*(1-sin(b(2)*x(:,1)*2*pi+b(3)))+b(4);
fit_stats=NaN(size(amp_list,2),4);
tval=linspace(0,time_max,samples);
tval=tval';

for idxa=1:size(amp_list,2)
    min_est=NaN(repeat_fit,2);
    amp_run=amp_list(idxa);
    for idxb=1:repeat_fit
        beta0=[amp_run,freq_truth,phase_truth,offset_truth];
        xval=modelfun(beta0,tval);
        xval=xval+normrnd(0,noise_amp,size(xval,1),1);
        names={'amp','freq','phase','offset'};

        beta0=beta0;
        mdl = fitnlm(tval,xval,modelfun,beta0,'CoefficientNames',names);
        min_est(idxb,:)=[mdl.Coefficients.Estimate(4),mdl.Coefficients.SE(4)];
        if fit_plots
            sfigure(1);
            set(gcf,'color','w')
            subplot(2,1,1)
            plot(tval,xval,'rx')
            hold on
            tpred=linspace(min(tval),max(tval),1e4)';
            xpred=predict(mdl,tpred);
            plot(tpred,xpred,'k')
            hold off
            subplot(2,1,2)
            plot(tval,predict(mdl,tval)-xval)
            pause(0.01)
        end       
end

fit_stats(idxa,:)=[std(min_est(:,1)),mean(min_est(:,2)),...
    mean(abs(min_est(:,1)-offset_truth)),rms(min_est(:,1)-offset_truth)];
fprintf('amp %.1f actual sd %.2e estimated %.2e mean error %.2e RMSE %.2e\n',...
    [amp_run,fit_stats(idxa,:)])

end
%%
figure(2)
clf;
raw_err=fit_stats(:,4)';
smoth_err=smoothdata(raw_err,'gaussian',10);
plot(amp_list,raw_err,'k')
hold on
plot(amp_list,smoth_err)
hold off
%%

%%
time_max=10;
tval=linspace(0,time_max,200);
tval=tval';
min_est=NaN(repeat_fit,2);
for idxb=1:10
    beta0=[amp_truth,dampt_truth,freq_truth,phase_truth,offset_truth,grad_truth];
    xval=modelfun(beta0,tval);
    xval=xval+normrnd(0,noise,size(xval,1),1);
    sfigure(1);
    plot(tval,xval,'rx')
    hold on
    names={'amp','dampt','freq','phase','offset','grad'};
    beta0=[amp_truth,dampt_truth,freq_truth,phase_truth,offset_truth,grad_truth];
    mdl = fitnlm(tval,xval,modelfun,beta0,'CoefficientNames',names);
    min_est(idxb,:)=[mdl.Coefficients.Estimate(3),mdl.Coefficients.SE(3)];
    tpred=linspace(min(tval),max(tval),1e4)';
    xpred=predict(mdl,tpred);
    plot(tpred,xpred,'k')
    hold off
end
fprintf('time %.1f actual sd %.2e estimated %.2e mean error %.2e\n',[time_max,std(min_est(:,1)),mean(min_est(:,2)),mean(abs(min_est(:,1)-freq_truth))])
fit_stats(idxa,:)=[std(min_est(:,1)),mean(min_est(:,2)),mean(abs(min_est(:,1)-freq_truth))];
