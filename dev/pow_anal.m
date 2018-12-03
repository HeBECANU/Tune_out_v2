function coefs = pow_anal(X,Y)
pow_fit_mdl = @(p,x) p(1) + p(2)*x.^p(3);
coefs = [0,6.5,-1.5];
% pow_fit = fitnlm(X,Y,pow_fit_mdl,coefs);
% coefs = pow_fit.Coefficients.Estimate;
plot(X,Y)
hold on
plot(X,pow_fit_mdl(coefs,X))
hold off

end