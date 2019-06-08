function out_st=global_fit_new(predictor,response,modelfun,gf_opt)

cost_fun=@(x) sqrt(nansum((modelfun(col_vec(x),predictor)-response).^2)/(numel(response)-numel(x)));

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
    gf_opt.level=3;
end
if ~isfield(gf_opt,'verbose')
    gf_opt.verbose=0;
end

out_st=[];
gf_opt.start=col_vec(gf_opt.start); %srarting parameters

lb=col_vec(min(gf_opt.domain,[],2));
ub=col_vec(max(gf_opt.domain,[],2));

if ~isequal(size(lb),size(ub),size(gf_opt.start))
    error('bounds and start are not the same size')
end

if sum(gf_opt.start>ub | gf_opt.start<lb)>0
    fprintf(2,'start is not in bounds\n') %red text
    if sum(gf_opt.start>ub)>0
        elm_idx=1:numel(gf_opt.start);
        elm_idx=elm_idx(gf_opt.start>ub);
        fprintf('indices exceeding upper bound %s\n',sprintf('%u,',elm_idx))
        fprintf('values %s\n',sprintf('%f,',gf_opt.start(elm_idx)))
    elseif sum(gf_opt.start<lb)>0
        elm_idx=1:numel(gf_opt.start);
        elm_idx=elm_idx(gf_opt.start<lb);
        fprintf('indices exceeding lower  bound %s\n',sprintf('%u,',elm_idx))
        fprintf('values %s\n',sprintf('%f,',gf_opt.start(elm_idx)))
    end
    error('start is not in bounds')
end


params=gf_opt.start;
cost_fun_inital=cost_fun(params);

out_st.rmse=inf;
out_st.params=params;
method_counter=0;
while out_st.rmse>gf_opt.rmse_thresh && method_counter<=gf_opt.level
    switch  method_counter
        case 0
            %initalize the output strucutre
            out_st.history.optimizer={'none'};
            out_st.history.cost=cost_fun_inital;
            out_st.history.best_params={params};
            out_st.params=params;
            out_st.rmse=cost_fun_inital;
        case 1
            out_st.history.optimizer{end+1}='fmincon';
            options = optimset('Display','off'); 
            params=fmincon(cost_fun,gf_opt.start,[],[],[],[],lb,ub,[],options);
            cost_fun_fmincon=cost_fun(params);
            params=col_vec(params);
            out_st.history.cost=[out_st.history.cost,cost_fun_fmincon];
            out_st.history.best_params{end+1}=params;
        case 2
           out_st.history.optimizer{end+1}='patternsearch';
            %use patternsearch to find the optima
            PSoptions = optimoptions(@patternsearch,... %'Display','iter'
            'MaxIterations',1e4,...
            'MaxFunctionEvaluations',1e5,...
            'MeshTolerance',1e-5,...
            'Display','off');
            [params,costend] = patternsearch(cost_fun,gf_opt.start,...
            [],[],[],[],...
            col_vec(min(gf_opt.domain,[],2)),col_vec(max(gf_opt.domain,[],2)),PSoptions);
            cost_fun_patternsearch=costend;
            params=col_vec(params);
            out_st.history.cost=[out_st.history.cost,cost_fun_patternsearch];
            out_st.history.best_params{end+1}=params;
        case 3
            out_st.history.optimizer{end+1}='godlike';
            % use godlike to find the optima
            [params,costend,outflag,details]=GODLIKE(cost_fun,...
                lb,ub,...
                [],...
                'AchieveFunVal',gf_opt.rmse_thresh*0.90,... %sometimes GODLIKE returns a value that is slightly higher
                'MaxFunEvals',1e5,...
                'MaxIters',100,...
                'display', 'off');
            params=col_vec(params);
            out_st.history.cost=[out_st.history.cost,costend];
            out_st.history.best_params{end+1}=params;
       case 4
            out_st.history.optimizer{end+1}='godlike_deep';
            % use godlike to find the optima
            [params,costend,outflag,details]=GODLIKE(cost_fun,...
                lb,ub,...
                [],...
                'AchieveFunVal',gf_opt.rmse_thresh*0.90,... %sometimes GODLIKE returns a value that is slightly higher
                'MaxFunEvals',1e6,...
                'MaxIters',400,...
                'display', 'off');
            params=col_vec(params);
            out_st.history.cost=[out_st.history.cost,costend];
            out_st.history.best_params{end+1}=params;     
    end
    if gf_opt.verbose>3
        fprintf('%s:%s out cost %f \n',mfilename,out_st.history.optimizer{end},out_st.history.cost(end))
    end
    method_counter=method_counter+1;
    [out_st.rmse,idx]=min(out_st.history.cost);
    out_st.params=out_st.history.best_params{idx};
    out_st.optimizer=out_st.history.optimizer{idx};
end

if gf_opt.verbose>3
        fprintf('%s:final out cost %f \n',mfilename,out_st.rmse)
end
    
if out_st.rmse>gf_opt.rmse_thresh
    warning('global_fit could not find a solution with rmse less than thresh(%f) min solution %f\n  ',gf_opt.rmse_thresh,out_st.rmse)
end



if  gf_opt.plot
    stfig('solution','add_stack',1);
    clf
    subplot(2,1,1)
    plot(predictor,response,'k')
    hold on
    plot(predictor,modelfun(gf_opt.start,predictor),'b')
    plot(predictor,modelfun(out_st.params,predictor),'r')
    hold off
    legend('data','inital guess','global fit')
    subplot(2,1,2)
    plot(predictor,modelfun(out_st.params,predictor)-response,'r')
    drawnow
end




end

