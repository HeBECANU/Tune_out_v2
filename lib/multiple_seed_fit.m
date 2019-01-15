function [optimum_params,optimum_ci,min_res] =  multiple_seed_fit(model,predictor,weights,target_val,intial_seed,max_seed,ub,lb)
%fitting procedure which uses multiple starting seeds
opt = optimset('MaxIter',1e5,...
            'UseParallel',0,'TolFun',1e-12,'TolX',1e-12,'MaxFunEvals',0.3e3,'display','off');
seed = intial_seed;
seed_num = 0;
min_res = 10;
p_num = length(seed);
id = 'MATLAB:singularMatrix';
warning('off',id)
while seed_num<max_seed && min_res>1e-5
    costfun = @(b) double(weights.*(model(b,predictor(:,1))-target_val))./double(sum(weights));
    [fitparam,resnorm,residual,exitflag,out_put,lm,J]=lsqnonlin(costfun,double(seed),double(lb),double(ub),opt);
    seed_num = seed_num + 1;
    %confidence interval
    ci = nlparci(fitparam,residual,'Jacobian',J); %estimates conffidence interval
    fiterror = ci(:,2)'-ci(:,1)';
    fiterror(isnan(fiterror))= 0.1;
    if resnorm<min_res
        optimum_params = fitparam;
        min_res = resnorm;
        optimum_ci = fiterror;
    end
    %select new seed
    if mod(seed_num,3) == 0
        seed =  optimum_params + (ub-lb).*0.15.*rand(1,p_num);
    elseif mod(seed_num,5) == 0
        seed = (ub-lb).*rand(1,p_num) + lb;
    else
        seed = fitparam + fiterror.*0.5.*(0.5-rand(1,p_num));
    end
    
    for ii = 1:p_num
        if seed(ii)>ub(ii)
            seed(ii) = ub(ii);
        elseif seed(ii)<lb(ii)
            seed(ii) = lb(ii);
        end
    end
end