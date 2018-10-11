%SIMULATE measuring atom number and outcoupling fraction by fitting
%compare fitting in semi log and fitting the exponential form directly
sim_det_num=zeros(300,1);
sim_rem_num=sim_det_num;
pulse_num=(1:size(sim_det_num,1))-1;
atoms_inital=1e6;
atoms=atoms_inital;
out_frac=5e-3;
QE=0.1;
for n=1:size(sim_det_num,1)
    sim_rem_num(n)=atoms;
    sim_det_num(n)=sum(rand(round(atoms*out_frac),1)<QE);
    atoms=atoms-atoms*out_frac;
end
%cant get the weights to work for the life of me
pred_det_noise=sqrt(sim_det_num);
pred_det_noise_log=pred_det_noise./sim_det_num;
weights=1./(pred_det_noise_log.^2);     
%,'Weights',weights
fit_mdl_log=fitlm(pulse_num,log(sim_det_num),'linear');  
xplot=linspace(min(pulse_num),max(pulse_num),1e3)';
[ymodel,ci]=predict(fit_mdl_log,xplot,'Prediction','observation','Alpha',0.03174);
plot(pulse_num,sim_det_num,'kx')
hold on
plot(xplot,exp(ymodel),'r-')
plot(pulse_num,sim_rem_num*QE*out_frac,'g-') %exp(ci(:,1))
plot(pulse_num,sim_rem_num*QE*out_frac-sqrt(sim_rem_num*QE*out_frac),'b-') %exp(ci(:,1))
plot(pulse_num,sim_rem_num*QE*out_frac+sqrt(sim_rem_num*QE*out_frac),'b-') %exp(ci(:,1))
plot(xplot,exp(ci(:,2)),'m-')
plot(xplot,exp(ci(:,1)),'m-')
hold off
set(gca, 'YScale', 'log')

fit1_res=[[-fit_mdl_log.Coefficients.Estimate(2),fit_mdl_log.Coefficients.SE(2)];...
    [-exp(fit_mdl_log.Coefficients.Estimate(1))/(fit_mdl_log.Coefficients.Estimate(2)*QE),...
    0]];
fit1_res(2,2)=(1/QE)*fit1_res(2,1)*sqrt(fit_mdl_log.Coefficients.SE(1)^2+...
    (fit_mdl_log.Coefficients.SE(2)/fit_mdl_log.Coefficients.Estimate(2))^2);
fprintf('semi log fit results\n')
fprintf('frac actual %.2e,predicted %.2e±%.2e\nerr sigma %.2f frac err act %.3f sd %.3f\n',...
    out_frac,fit1_res(1,1),fit1_res(1,2),(fit1_res(1,1)-out_frac)/fit1_res(1,2),...
    (fit1_res(1,1)-out_frac)/out_frac,fit1_res(1,2)/out_frac)
fprintf('no actual %.2e,predicted %.2e±%.2e\nerr sigma %.2f frac err act %.3f sd %.3f\n',...
    atoms_inital,fit1_res(2,1),fit1_res(2,2),(fit1_res(2,1)-atoms_inital)/fit1_res(2,2)...
    ,(fit1_res(2,1)-atoms_inital)/atoms_inital,fit1_res(2,2)/atoms_inital)

%lets compare that to doing a fit of a nonlinear model

%%
fo = statset('TolFun',1e-30,...
    'TolX',1e-50,...
    'MaxIter',1e7,...
    'UseParallel',1); %'DerivStep',1e-5
modelfun = @(b,x) b(1)*b(2)*((1-b(2)).^(x(:,1)-1));
beta0=[1e6,1e-2];
pred_det_noise=sqrt(sim_det_num);
weights=1./(pred_det_noise.^2);
%,'Weights',weights 
fit_mdl_pow=fitnlm(pulse_num,sim_det_num,modelfun,...
   beta0,'CoefficientNames',{'N0','eta'},'Options',fo);  
xplot=linspace(min(pulse_num),max(pulse_num),1e3)';
[ymodel,ci]=predict(fit_mdl_pow,xplot,'Prediction','observation','Alpha',0.3174);
plot(pulse_num,sim_det_num,'kx')
hold on
plot(xplot,ymodel,'r-')
plot(pulse_num,sim_rem_num*QE*out_frac,'g-') %exp(ci(:,1))
plot(pulse_num,sim_rem_num*QE*out_frac-sqrt(sim_rem_num*QE*out_frac),'b-') %exp(ci(:,1))
plot(pulse_num,sim_rem_num*QE*out_frac+sqrt(sim_rem_num*QE*out_frac),'b-') %exp(ci(:,1))
plot(xplot,ci(:,2),'m-')
plot(xplot,ci(:,1),'m-')
hold off
set(gca, 'YScale', 'log')

fit2_res=[[fit_mdl_pow.Coefficients.Estimate(2),fit_mdl_pow.Coefficients.SE(2)];...
    [fit_mdl_pow.Coefficients.Estimate(1)/QE,...
    fit_mdl_pow.Coefficients.SE(1)/QE]];
fprintf('linear power fit results\n')
fprintf('frac actual %.2e,predicted %.2e±%.2e\nerr sigma %.2f frac err act %.3f sd %.3f\n',...
    out_frac,fit2_res(1,1),fit2_res(1,2),(fit2_res(1,1)-out_frac)/fit2_res(1,2),...
    (fit2_res(1,1)-out_frac)/out_frac,fit2_res(1,2)/out_frac)
fprintf('no actual %.2e,predicted %.2e±%.2e\nerr sigma %.2f frac err act %.3f sd %.3f\n',...
    atoms_inital,fit2_res(2,1),fit2_res(2,2),(fit2_res(2,1)-atoms_inital)/fit2_res(2,2)...
    ,(fit2_res(2,1)-atoms_inital)/atoms_inital,fit2_res(2,2)/atoms_inital)



%data_14param{1}{5}{1}.num_in_win