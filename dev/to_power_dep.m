 to_vals_lin_quad = [
     0.8,725736332.4,127,725736120.4,157;
     0.2,725736044.9,182,725736234.7,224;
     0.7,725735895.5,100,725735679.2,58;
     0.4,725736096.6,154,725736029.7,153;
     0.6,725735898.6,159,725736124.8,156];
 
 
 %fit and plot

 mean_to_freq=mean(to_vals_lin_quad(:,2));
 errorbar(to_vals_lin_quad(:,1),to_vals_lin_quad(:,2)-mean_to_freq,to_vals_lin_quad(:,3),'x')
 xdat=to_vals_lin_quad(:,1);
 ydat=to_vals_lin_quad(:,2)-mean_to_freq;
 yunc=to_vals_lin_quad(:,3);

     
 
modelfun = @(b,x) b(1)+x(:,1).*b(2);
opts = statset('nlinfit');
beta0 = [0,1]; %intial guesses 
%'Weights',wlin,
fit_mdl_lin = fitnlm(xdat,ydat,modelfun,beta0,...
    'Options',opts,'CoefficientNames' ,{'offset','lin'})

paddx=0.1;
xmdl_samp=linspace(min(xdat)-paddx,max(xdat)+paddx,1e3)';
ci_plot=1-erf(1/sqrt(2));%1-erf(zvalue/sqrt(2)) %confidence interval for cutting outliers
[ymdl_val,ymdl_ci]=predict(fit_mdl_lin,xmdl_samp,'Prediction','observation','Alpha',ci_plot);

sfigure(701)
errorbar(to_vals_lin_quad(:,1),ydat,yunc,'x')
hold on
plot(xmdl_samp,ymdl_val,'b-','LineWidth',1.6)
plot(xmdl_samp,ymdl_ci(:,1),'r-','LineWidth',1.6)
plot(xmdl_samp,ymdl_ci(:,2),'r-','LineWidth',1.6)
hold off