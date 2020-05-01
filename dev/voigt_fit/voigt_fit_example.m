% voigt fit play


% addpath('../../lib/Core_BEC_Analysis/lib/') %add the path to set_up_project_path
% set_up_project_path('../..')



%% lets just check that we can make some voit pofiles
% xcross check a direct convolution method
% against 2 calculation methods in voigt_function_1d

w_gauss=0.1;
w_lorz=0.2;
center=0.0;

xsamp_comp=linspace(-10,10,1e5+1); %odd number of bins allows better results from conv
xsamp_comp=col_vec(xsamp_comp);
dx_comp=abs(xsamp_comp(2)-xsamp_comp(1));

gauss_component=gaussian_function_1d(xsamp_comp,w_gauss,center,'norm','int');
trapz(dx_comp,gauss_component)
lorentzian_component=lorentzian_function_1d(xsamp_comp,w_lorz,center,'norm','int');
trapz(dx_comp,lorentzian_component)

stfig('plot componets');
clf
subplot(2,1,1)
plot(xsamp_comp,gauss_component)
hold on
plot(xsamp_comp,lorentzian_component)

conv_comp=conv(gauss_component,lorentzian_component,'same')*dx_comp;

xsamp_conv=xsamp_comp(1)+dx_comp*(0:(numel(conv_comp)-1));
xsamp_conv=col_vec(xsamp_conv)

plot(xsamp_conv,conv_comp)

tic
vo_approx=voigt_function_1d(xsamp_conv,w_gauss,w_lorz,center,'method','approx');
toc
plot(xsamp_conv,vo_approx,'--')


vo_fadeva=voigt_function_1d(xsamp_conv,w_gauss,w_lorz,center,'method','fadd');

plot(xsamp_conv,vo_fadeva,'-.')

legend('gauss','lorz','conv','approx',' fadeva')
hold off

approx_vfwhm=voigt_approx_fwhm(w_gauss,w_lorz);

xlim([-1,1]*1*approx_vfwhm+center)


subplot(2,1,2)
error_approx=conv_comp-vo_approx;
plot(xsamp_conv,error_approx)
hold on
error_fadeva=conv_comp-vo_fadeva;
plot(xsamp_conv,error_fadeva,'--')
hold off
xlim([-1,1]*1*approx_vfwhm+center)

% that all seems ok


%%
shaded_curve_ci_lines=false;
fit_fun = @(b,x) voigt_function_1d(x,b(1),b(2),b(3),b(4),b(5),'norm','amp','method','approx');
%fit_fun = @(b,x) gaussian_function_1d(xsamp_conv,b(1),b(3),b(4),b(5),'norm','amp');
% this does not work with the predict function
%fit_fun1 = @(b,x) test_fun_1(x,b(1),b(3),b(4),b(5))+b(2)*0;

%fit_fun=@(b,x) b(4)*exp(-(1/2)*((x-b(3))./b(1)).^2)+b(5)+b(2)*0;


cof_names={'sigma','gamma','mu','amp','offset'};
opt = statset('TolFun',1e-10,'TolX',1e-10,...
    'MaxIter',1e4,... %1e4
    'UseParallel',1);
beta0=[0.5,0.5,0.5,1.1,0];

predictor=xsamp_conv;
response=vo_fadeva+randn(size(predictor))*0.1;
% weights=ones(size(predictor));
% %'Weights',weights

fitobj=fitnlm(predictor,response,fit_fun,beta0,...
    'options',opt,...
    'CoefficientNames',cof_names)

%%
stfig('fit results')
clf
plot(predictor,response,'k.')
hold on


tplotvalues=linspace(min(predictor),max(predictor),1e4);
tplotvalues=col_vec(tplotvalues);
%[prediction,ci]=predict(fitobject,tplotvalues); %'Alpha',1-erf(1/sqrt(2)),'Prediction','observation'
[amp_pred,ci_obs]=fitobj.predict(tplotvalues,'Alpha',1-erf(1/sqrt(2)),'Prediction','observation'); %
[~,ci_curve]=fitobj.predict(tplotvalues,'Alpha',1-erf(1/sqrt(2)),'Prediction','curve'); %
size(amp_pred)
size(tplotvalues)
%
shaded_curve_ci_lines=false;
if shaded_curve_ci_lines
    color_shaded=[0.9,1,1];
    patch([tplotvalues', fliplr(tplotvalues')], [ci_obs(:,1)', fliplr(ci_obs(:,2)')], color_shaded,'EdgeColor','none')  %[1,1,1]*0.80
    hold on
else
    color_shaded=[0,1,0];
    plot(tplotvalues,ci_obs(:,1),'-','LineWidth',1.5,'Color',color_shaded)
    hold on
    plot(tplotvalues,ci_obs(:,2),'-','LineWidth',1.5,'Color',color_shaded)
end  

shaded_curve_ci_lines=false;
if shaded_curve_ci_lines
    color_shaded=[0.9,1,1];
    patch([tplotvalues', fliplr(tplotvalues')], [ci_curve(:,1)', fliplr(ci_curve(:,2)')], color_shaded,'EdgeColor','none')  %[1,1,1]*0.80
    hold on
else
    color_shaded=[0,1,0];
    plot(tplotvalues,ci_curve(:,1),'-','LineWidth',1.5,'Color',color_shaded)
    hold on
    plot(tplotvalues,ci_curve(:,2),'-','LineWidth',1.5,'Color',color_shaded)
end  

plot(tplotvalues,amp_pred,'r-','LineWidth',1.0)




hold off 


