function out_st=global_fit(predictor,response,modelfun,rf_opt)

cost_fun=@(x) sqrt(sum((modelfun(col_vec(x),predictor)-response).^2)/(numel(response)-numel(x)));

%inital start rmse
%cost_fun(rf_opt.start)
%cost_fun(rf_opt.start+1)
% plot(predictor,modelfun(rf_opt.start,predictor))
% hold on
% plot(predictor,response)

if ~isfield(rf_opt,'plot')
    rf_opt.plot=false;
end


out_st=[];

lb=col_vec(min(rf_opt.domain,[],2));
ub=col_vec(max(rf_opt.domain,[],2));
params=rf_opt.start;
if cost_fun(params)>rf_opt.rmse_thresh
    options = optimset('Display','off'); 
    params=fmincon(cost_fun,rf_opt.start,[],[],[],[],lb,ub,[],options);
    if cost_fun(params)>rf_opt.rmse_thresh
        % use patternsearch to find the optima
        PSoptions = optimoptions(@patternsearch,... %'Display','iter'
            'MaxIterations',1e4,...
            'MaxFunctionEvaluations',1e5,...
            'MeshTolerance',1e-5,...
            'Display','off');
        [params,costend] = patternsearch(cost_fun,rf_opt.start,...
            [],[],[],[],...
            col_vec(min(rf_opt.domain,[],2)),col_vec(max(rf_opt.domain,[],2)),PSoptions);

        if costend>rf_opt.rmse_thresh
            % use godlike to find the optima
            [params,costend,outflag,details]=GODLIKE(cost_fun,...
                lb,ub,...
                [],...
                'AchieveFunVal',rf_opt.rmse_thresh*0.90,... %sometimes GODLIKE returns a value that is slightly higher
                'MaxFunEvals',1e5,...
                'MaxIters',100,...
                'display', 'off');
            if costend>rf_opt.rmse_thresh
                warning('global_fit could not find a solution with rmse less than thresh(%f) min solution %f\n    ',rf_opt.rmse_thresh,costend);
            end
        end
    end
end


if  rf_opt.plot
    stfig('solution','add_stack',1);
    clf
    subplot(2,1,1)
    plot(predictor,response,'k')
    hold on
    plot(predictor,modelfun(params,predictor),'r')
    subplot(2,1,2)
    plot(predictor,modelfun(params,predictor)-response,'r')
    drawnow
end

out_st.params=params;
out_st.rmse=cost_fun(params);



end

