%test_fit_poly_with_int
% generate some fake data and try to fit to it


grad=rand(1)*10;
curve=rand(1);
yintercept=rand(1)*10;
points=100;
noise_level_y=2;
noise_level_x=0;

%dep_fun=@(x) yintercept+grad.*x+curve.*x.^2;
dep_fun=@(x) yintercept+grad.*x;
xdat_no_noise=rand(points,1)*10;
ydat_no_noise=dep_fun(xdat_no_noise);



%%
iimax=100;
fit_order=2;

fit_xint=[];
fit_xint_unc=[];
do_plot=false;
for ii=1:iimax
    ydat=ydat_no_noise+randn(numel(xdat_no_noise),1)*noise_level_y;
    xdat=xdat_no_noise+randn(numel(xdat_no_noise),1)*noise_level_x;
    wdat=xdat_no_noise*0+1;
    xywdat=num2cell([xdat,ydat,wdat],2);
    out=fit_poly_data(xywdat,fit_order);
    fit_xint(ii)=out.x_intercept.val;
    fit_xint_unc(ii)=out.x_intercept.unc.with_cov;
    if do_plot
        stfig('eval fit performance')
        clf
        plot(xdat_no_noise,ydat,'x')
        xsamp=col_vec(linspace(min(xdat_no_noise),max(xdat_no_noise),1e4));
        [yfit,yci]=predict(fit_mdl{ii},xsamp,'Prediction','observation','Alpha',1-erf(1/sqrt(2)));
        hold on
        plot(xsamp,yfit)
        plot(xsamp,yci(:,1))
        plot(xsamp,yci(:,2))
        hold off
        drawnow
    end
end


hist(fit_xint)
mean_xfit_int=mean(fit_xint)
acutal_xint=fzero(dep_fun,0)
fprintf('mean x intercept %f, true value %f\n',mean_xfit_int,acutal_xint)
fprintf('std x intercept %f, mean fit unc %f\n',std(fit_xint),mean(fit_xint_unc))
fprintf('mean-true in std %f \n',(mean_xfit_int-acutal_xint)/std(fit_xint))


