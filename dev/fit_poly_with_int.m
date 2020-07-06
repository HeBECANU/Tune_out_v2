function out=fit_poly_with_int(in,poly_order,use_robust,use_weights)
% this function shoul fit for the x intercept of a polynomial regression
% TODO
% is it better to turn fit on its side to get intercept
% test the polyval function
% add in the coefficinet names
% decent inital guesses

% see fit_poly_with_int
if iscell(in)
    inmat=cell2mat(in);
else
    inmat=in;
end
xdat=inmat(:,1);
ydat=inmat(:,2);
wdat=inmat(:,3);

modelfun=@(b,x) polyval(fliplr(b),x);
%modelfun = @(b,x) (repmat(x(:,1),1,n+1).^(0:n))*(b(1:(n+1))'); %simple linear model
%
opts = statset('nlinfit');
%robust fit
if use_robust
    opts.RobustWgtFun = 'welsch' ; %a bit of robust fitting
    opts.Tune = 1;
end
% if robust cant use the combined err model
%,'ErrorModel','combined'
beta0 = [1e-5,1e-2.*ones(1,poly_order)]; %intial guesses
if use_weights
    fit_mdl = fitnlm(xdat,ydat,modelfun,beta0,'Options',opts,'Weight',wdat);
else
    fit_mdl = fitnlm(xdat,ydat,modelfun,beta0,'Options',opts);
end


% to find the x intercept we find the roots of the polynomial
%fzero(@(x) predict(fit_mdl,x),0)
mdl_zeros = roots(fliplr(fit_mdl.Coefficients.Estimate(:)'));
mdl_zeros=mdl_zeros(imag(mdl_zeros)==0);
if numel(mdl_zeros)<1
    warning('fit has found no roots of polynomial expression')
    x_intercept=nan;
else
    linear_xintercept=-fit_mdl.Coefficients.Estimate(1)/fit_mdl.Coefficients.Estimate(2);
    [~, indx] = min(abs(mdl_zeros-linear_xintercept));
    x_intercept = mdl_zeros(indx);
end
mdl_coef=fit_mdl.Coefficients;
covm=fit_mdl.CoefficientCovariance;
 
if poly_order==1 %propogation for the linear intercept
    unc_x_int_no_cov=abs((mdl_coef.Estimate(1)/mdl_coef.Estimate(2)))*...
        sqrt(...
        (mdl_coef.SE(1)/mdl_coef.Estimate(1))^2+...
        (mdl_coef.SE(2)/mdl_coef.Estimate(2))^2 ...
        );

    %calculate the error propagation with the covariance included
    %generaly makes a tiny change ~1e-3
    unc_x_int_cov=abs((mdl_coef.Estimate(1)/mdl_coef.Estimate(2)))*...
        sqrt(...
        (mdl_coef.SE(1)/mdl_coef.Estimate(1))^2+...
        (mdl_coef.SE(2)/mdl_coef.Estimate(2))^2 ...
        - 2*covm(2,1)/(mdl_coef.Estimate(1)*mdl_coef.Estimate(2))...
        );
elseif poly_order==2 %propogation for quadratic intercept
    det = sqrt(mdl_coef.Estimate(2)^2-4*mdl_coef.Estimate(1)*mdl_coef.Estimate(3));
    c = mdl_coef.Estimate(1);
    b = mdl_coef.Estimate(2);
    a= mdl_coef.Estimate(3);
    unc_a_c_2 = (mdl_coef.SE(1)^2+(2*c*a+b*(det-b))^2/(4*a^4)*mdl_coef.SE(3)^2)*1/det^2;
    unc_b_2 = (b/det-1)^2*1/(4*a^2)*mdl_coef.SE(2)^2;
    unc_x_int_no_cov=sqrt(...
        unc_a_c_2+...
        unc_b_2 );

    %calculate the error propagation with the covariance included
    %generaly makes a tiny change ~1e-3
    unc_ab = (-b^3+3*a*b*c+det*(b^2-a*c))*covm(3,2);
    unc_ac = (a*b^2-2*a^2*c-a*b*det)*covm(3,1);
    unc_bc = (-a^2*b+a^2*det)*covm(2,1);
    unc_x_int_cov=sqrt(unc_x_int_no_cov^2+...
    1/(a^3*det^2)*(unc_ab+unc_ac+unc_bc)...
                        );
else
    unc_x_int_no_cov=nan;
    unc_x_int_cov=nan;
    warning('x intercept error propagation has not yet been implemented for this order')
end

out.x_intercept.val=x_intercept;
out.x_intercept.unc.with_cov=unc_x_int_cov;
out.x_intercept.unc.no_cov=unc_x_int_no_cov;
out.fit_mdl=fit_mdl;


end



% function out=poly_mdl(b,x)
% %a=1:2;
% % x_pow=(repmat(x(:,1),1,n+1).^(0:n));
% % out= x_pow*(b(1:(n+1))');
% 
% end