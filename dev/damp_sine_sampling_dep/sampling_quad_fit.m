%sampling errors linfit
%studying if there is an optimal rang to sample a inear fucntion (with noise in y) in oder to determine the x axis
%crossing

% BEGIN USER VAR-------------------------------------------------
grad_truth=1;
xintercept_truth=1e-3; %x value of the x axis intercept
quad_truth=5e-5;
noise_amp=1;
repeat_fit=1e1;
range_list=logspace(0,2,1e2);
samples=300;
plot_fits=false;
% END USER VAR-----------------------------------------------------------

quadmodel=@(b,x) b(1)+b(2)*x+b(3)*x.^2;
yintercept_truth=-xintercept_truth*grad_truth; %transform into y intercept value
est_fit_param=nan(repeat_fit,3);
denom_se_grad=[];
iimax=size(range_list,2);
fit_stats=nan(iimax,size(est_fit_param,2),3);
beta0=[yintercept_truth,grad_truth,quad_truth];

fprintf('%04i:%04i',iimax,0)
for idxa=1:iimax
    fprintf('\b\b\b\b%04i',idxa)
    range=range_list(idxa);
    xval=linspace(-range,range,samples)';
    for idxb=1:repeat_fit
        ydat=quadmodel(beta0,xval);
        ydat=ydat+normrnd(0,noise_amp,size(ydat,1),1);
        opts = statset('TolFun',1e-10,'TolX',1e-10); %,'fitnlm'
        fit_mdl = fitnlm(xval,ydat,quadmodel,beta0,'Options',opts);
        %xint_samp=-fit_mdl.Coefficients.Estimate(2)/fit_mdl.Coefficients.Estimate(1);
        xint_samp= fzero(@(x) quadmodel(fit_mdl.Coefficients.Estimate,x),0);
        est_fit_param(idxb,:)=[fit_mdl.Coefficients.Estimate(2),xint_samp,fit_mdl.Coefficients.Estimate(1)];
        if plot_fits
            sfigure(1);
            set(gcf,'color','w')
            subplot(2,1,1)
            plot(xval,ydat,'rx')
            hold on
            xplot=linspace(min(xval),max(xval),1e4)';
            plot(xplot,predict(fit_mdl,xplot),'k')
            plot(xval,quadmodel(beta0,xval),'g-')
            hold off
            legend('sampled data','fitted model','ground truth model')
            subplot(2,1,2)
            plot(xval,predict(fit_mdl,xval)-ydat,'x')
            ylabel('residuals')
            pause(0.01)
        end       
    end
    fit_stats(idxa,:,:)=[mean(est_fit_param)',std(est_fit_param)',median(est_fit_param)'];
    denom_se_grad(idxa)=sqrt(sum((xval-mean(xval)).^2));
end
fprintf('\n')


%%
%https://stats.stackexchange.com/questions/289457/proof-for-the-standard-error-of-parameters-in-linear-regression
predicted_yintercept_unc=noise_amp/sqrt(samples);
%predicted_xintercept_unc=abs(fit_stats(:,2,1)).*sqrt((fit_stats(:,3,2)./fit_stats(:,3,1)).^2+(fit_stats(:,1,2)./fit_stats(:,1,1)).^2);
predicted_grad_unc=noise_amp./denom_se_grad';
predicted_xintercept_unc=abs(xintercept_truth).*sqrt((predicted_yintercept_unc./yintercept_truth).^2 ...
                            +(predicted_grad_unc./grad_truth).^2);

%yintercept value , or fit_stats(:,3,1)
% value abs(fit_stats(:,2,1)) or  xintercept_truth
figure(2)
clf
set(gcf,'color','w')
subplot(3,3,1)
semilogx(range_list,fit_stats(:,1,1))
xl=xlim;
line(xl,[1,1]*grad_truth,'Color','k')
title('Gradient')
ylabel('mean gradient')
xlabel('range')
subplot(3,3,2)
semilogx(range_list,fit_stats(:,3,1))
xl=xlim;
line(xl,[1,1]*yintercept_truth,'Color','k')
ylabel('mean y intercept')
title('Y intercept')
xlabel('range')
subplot(3,3,3)
semilogx(range_list,fit_stats(:,2,1))
xl=xlim;
line(xl,[1,1]*xintercept_truth,'Color','k')
ylabel('mean x intercept')
title('X intercept')
xlabel('range')
subplot(3,3,4)
semilogx(range_list,fit_stats(:,1,3))
xl=xlim;
line(xl,[1,1]*grad_truth,'Color','k')
ylabel('median gradient')
xlabel('range')
subplot(3,3,5)
semilogx(range_list,fit_stats(:,3,3))
xl=xlim;
line(xl,[1,1]*yintercept_truth,'Color','k')
ylabel('median y intercept')
xlabel('range')
subplot(3,3,6)
semilogx(range_list,fit_stats(:,2,3))
xl=xlim;
line(xl,[1,1]*xintercept_truth,'Color','k')
ylabel('median x intercept')
xlabel('range')
subplot(3,3,7)
semilogx(range_list,fit_stats(:,1,2))
hold on
semilogx(range_list,predicted_grad_unc,'k')
hold off
ylabel('std gradient')
xlabel('range')
subplot(3,3,8)
semilogx(range_list,fit_stats(:,3,2))
xl=xlim;
line(xl,[1,1]*predicted_yintercept_unc,'Color','k')
ylabel('std y intercept')
xlabel('range')
subplot(3,3,9)
semilogx(range_list,fit_stats(:,2,2))
hold on
semilogx(range_list,predicted_xintercept_unc,'k')
hold off
ylabel('std x intercept')
xlabel('range')

