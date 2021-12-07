%try and find the optimal time to sample a sine with some decay
%constant

amp_truth   =2.8;
dampt_truth =0.665;   
freq_truth =420.4;
phase_truth =0;
offset_truth=0;
grad_truth  =0;
noise=0.64;
repeat_fit=3e1;
time_max_list=linspace(0.4,4,5e2);
freq_guess=45.4;
num=1e5;
samp_time=8e-3;
noise_no_avg=0.64*sqrt(num/300);

modelfun=@(b,x) b(1)*exp(-x(:,1)/max(b(2),0.01)).*sin(b(3)*x(:,1)*2*pi+b(4))+b(5)+b(6)*x(:,1);
fit_stats=NaN(size(time_max_list,2),4);

fit_plots=false;
for idxa=1:size(time_max_list,2)
    time_max=time_max_list(idxa);
    pulses=round(time_max/samp_time);
    tval=(1:(pulses-1))*samp_time;
    tval=tval';
    freq_est=NaN(repeat_fit,3);
    parfor idxb=1:repeat_fit
        beta0=[amp_truth,dampt_truth,freq_truth,phase_truth,offset_truth,grad_truth];
        xval=modelfun(beta0,tval);
        xval=xval+normrnd(0,noise_no_avg/sqrt(num/pulses),size(xval,1),1);
        names={'amp','dampt','freq','phase','offset','grad'};

        beta0=[amp_truth,dampt_truth,freq_guess,phase_truth,offset_truth,grad_truth];
        mdl = fitnlm(tval,xval,modelfun,beta0,'CoefficientNames',names);
        freq_est(idxb,:)=[mdl.Coefficients.Estimate(3),mdl.Coefficients.SE(3),3/samp_time+mdl.Coefficients.Estimate(3)];
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

fit_stats(idxa,:)=[std(freq_est(:,3)),mean(freq_est(:,2)),...
    mean(abs(freq_est(:,3)-freq_truth)),rms(freq_est(:,3)-freq_truth)];
fprintf('time %.1f actual sd %.2e estimated %.2e mean error %.2e RMSE %.2e\n',...
    [time_max,fit_stats(idxa,:)])

end
figure(2)
clf;
smoth_err=smoothdata(fit_stats(:,4),'gaussian',20);
plot(time_max_list,smoth_err)

%%
figure(2)
clf;
smoth_err=smoothdata(fit_stats(:,4),'gaussian',40);
plot(time_max_list,smoth_err)


%%
time_max=10;
tval=linspace(0,time_max,200);
tval=tval';
freq_est=NaN(repeat_fit,2);
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
    freq_est(idxb,:)=[mdl.Coefficients.Estimate(3),mdl.Coefficients.SE(3)];
    tpred=linspace(min(tval),max(tval),1e4)';
    xpred=predict(mdl,tpred);
    plot(tpred,xpred,'k')
    hold off
end
fprintf('time %.1f actual sd %.2e estimated %.2e mean error %.2e\n',[time_max,std(freq_est(:,1)),mean(freq_est(:,2)),mean(abs(freq_est(:,1)-freq_truth))])
fit_stats(idxa,:)=[std(freq_est(:,1)),mean(freq_est(:,2)),mean(abs(freq_est(:,1)-freq_truth))];
