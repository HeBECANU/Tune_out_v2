function out_st=global_fit(predictor,response,modelfun,gf_opt)

cost_fun=@(x) sqrt(sum((modelfun(col_vec(x),predictor)-response).^2)/(numel(response)-numel(x)));

%inital start rmse
%cost_fun(rf_opt.start)
%cost_fun(rf_opt.start+1)
% plot(predictor,modelfun(rf_opt.start,predictor))
% hold on
% plot(predictor,response)

if ~isfield(gf_opt,'plot')
    gf_opt.plot=false;
end
if ~isfield(gf_opt,'level')
    gf_opt.level=5;
end


out_st=[];

lb=col_vec(min(gf_opt.domain,[],2));
ub=col_vec(max(gf_opt.domain,[],2));
params=gf_opt.start;
out_st.optimizer='none';
if cost_fun(params)>gf_opt.rmse_thresh && gf_opt.level>0
    options = optimset('Display','off'); 
    params=fmincon(cost_fun,gf_opt.start,[],[],[],[],lb,ub,[],options);
    out_st.optimizer='fmincon';
    if cost_fun(params)>gf_opt.rmse_thresh && gf_opt.level>1
        % use patternsearch to find the optima
        PSoptions = optimoptions(@patternsearch,... %'Display','iter'
            'MaxIterations',1e4,...
            'MaxFunctionEvaluations',1e5,...
            'MeshTolerance',1e-5,...
            'Display','off');
        [params,costend] = patternsearch(cost_fun,gf_opt.start,...
            [],[],[],[],...
            col_vec(min(gf_opt.domain,[],2)),col_vec(max(gf_opt.domain,[],2)),PSoptions);
        out_st.optimizer='patternsearch';
        if costend>gf_opt.rmse_thresh && gf_opt.level>2
            % use godlike to find the optima
            [params,costend,outflag,details]=GODLIKE(cost_fun,...
                lb,ub,...
                [],...
                'AchieveFunVal',gf_opt.rmse_thresh*0.90,... %sometimes GODLIKE returns a value that is slightly higher
                'MaxFunEvals',1e5,...
                'MaxIters',100,...
                'display', 'off');
            out_st.optimizer='GODLIKE';
            if costend>gf_opt.rmse_thresh
                warning('global_fit could not find a solution with rmse less than thresh(%f) min solution %f\n    ',gf_opt.rmse_thresh,costend);
            end
        end
    end
end


if  gf_opt.plot
    stfig('solution','add_stack',1);
    clf
    subplot(2,1,1)
    plot(predictor,response,'k')
    hold on
    plot(predictor,modelfun(gf_opt.start,predictor),'b')
    plot(predictor,modelfun(params,predictor),'r')
    hold off
    legend('data','inital guess','global fit')
    subplot(2,1,2)
    plot(predictor,modelfun(params,predictor)-response,'r')
    drawnow
end

out_st.params=params;
out_st.rmse=cost_fun(params);



end

