xdata = [280,340,0,50,140,210];
%ydata = [20.9,60.9,141,270,33.2,293];
ydata = [59.6,107,130,126,235,160];
beta0 = [100,1,max(ydata)/2];
fitnlm(xdata,ydata,mdl,beta0)
mdl = @(b,x) b(1).*sin(2.*x(:,1).*pi/180+b(2))+b(3);
fit_mdl = fitnlm(xdata,ydata,mdl,beta0);
sfigure(23928374);
clf
X=linspace(0,360,100);plot(X,mdl(fit_mdl.Coefficients.Estimate,X'))
hold on
scatter(xdata,ydata,'kx')
xlabel('pol angle')
ylabel('power transmitted')