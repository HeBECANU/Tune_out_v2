function out_polz_state=pol_data_query_basic(pol_opts)
% a partial rewrite of pol_data_import
%% the idea of this function is to create a bunch of differenet models of the pol state for hwp/ qwp data angles
% Inputs
% pol_opts.location = 'pre_right';%post, pre_cen, pre_left, pre_right
% pol_opts.predict_method = 'fit';%'interp'; %obs (obsovation) fit (pertial fit) full_fit (fit with all parameters free)
% where you want to query a given polarization model 
% pol_opts.hwp =  use nan to mean no waveplate
% pol_opts.qwp = 

if ~isequal(size(pol_opts.hwp),size(pol_opts.qwp))
    error('query hwp and qwp not the same size')
end

out_polz_state=[];
out_polz_state.v.val=pol_opts.hwp*nan;
out_polz_state.v.unc=pol_opts.hwp*nan;
out_polz_state.theta.val=pol_opts.hwp*nan;
out_polz_state.theta.unc=pol_opts.hwp*nan;
out_polz_state.cont.val=pol_opts.hwp*nan;
out_polz_state.cont.unc=pol_opts.hwp*nan;

sin_mdl = @(b,x) b(1).*sin(b(2).*x.*pi/180+b(3))+b(4);
sin_mdl_fixed_freq4 = @(b,x) b(1).*sin(4.*x.*pi/180+b(2))+b(3);
sin_mdl_fixed_freq2 = @(b,x) b(1).*sin(2.*x.*pi/180+b(2))+b(3);
sin_mdl_abs = @(b,x) abs(b(1).*sin(b(2).*x.*pi/180+b(3))+b(4));
lin_mdl = @(b,x) b(1).*x+b(2);

post_rotation_angle=90;
pre_rotation_angle=-90;
pre_angle_flip=-1;
pre_hand_flip=1;


if strcmp(pol_opts.location,'post')
    pol_data_table=readtable('.\data\polz_data_post_window.csv');
    %pol_data_val.Properties.VariableNames
    %'qwp_angle_deg,hwp_angle_deg,max_power_uw,max_power_angle_deg,min_power_uw,min_power_angle_deg,handedness,'
    %if we want to use the observation method
    
    %to keep things in a consistent reference frame rotate the angles so that they are relative to up
    pol_data_table.min_power_angle_deg=pol_data_table.min_power_angle_deg+post_rotation_angle;
    
    pol_v =pol_data_table.handedness.*2.*sqrt(pol_data_table.max_power_uw.*pol_data_table.min_power_uw)...
        ./(pol_data_table.max_power_uw+pol_data_table.min_power_uw);%the V parameter for each run
    pol_cont = (pol_data_table.max_power_uw-pol_data_table.min_power_uw)...
        ./(pol_data_table.max_power_uw+pol_data_table.min_power_uw);%contrast
    pol_theta = pol_data_table.min_power_angle_deg.*pi/180; %using min pow angle
    
    
elseif strcmp(pol_opts.location,'pre_cen')
    % data 2019-03-20 in lab book
    pol_data_table_pre=readtable('.\data\polz_data_pre_window_cen.csv');
    %pol_data_val.Properties.VariableNames
    %'qwp_angle_deg,hwp_angle_deg,max_power_uw,max_power_angle_deg,min_power_uw,min_power_angle_deg,handedness,'
    pol_data_table_pre.min_power_angle_deg=pre_angle_flip*(pol_data_table_pre.min_power_angle_deg+pre_rotation_angle);
    pol_data_table_pre.handedness=pol_data_table_pre.handedness*pre_hand_flip;
    mean_power=mean(pol_data_table_pre.max_power_uw+pol_data_table_pre.min_power_uw);
    fprintf('mean power in pre_cen %f\n',mean_power)
    
    
    pol_data_table_post=readtable('.\data\polz_data_post_window.csv');
    pol_data_table_post.min_power_angle_deg=pol_data_table_post.min_power_angle_deg+post_rotation_angle;
    %to keep things in a consistent reference frame rotate the angles so that they are relative to up
    %todo, do for all the other angles
    % use the qwp data from this
    mask_post_data_to_use=~isnan(pol_data_table_post.qwp_angle_deg);
    pol_data_table_post=pol_data_table_post(mask_post_data_to_use,:);
    
    commom_var_names=intersect(pol_data_table_post.Properties.VariableNames,...
                                pol_data_table_pre.Properties.VariableNames);
    
    pol_data_table=[pol_data_table_pre(:,commom_var_names);pol_data_table_post(:,commom_var_names)];
    
    
    %if we want to use the observation method
    pol_v =pol_data_table.handedness.*2.*sqrt(pol_data_table.max_power_uw.*pol_data_table.min_power_uw)...
        ./(pol_data_table.max_power_uw+pol_data_table.min_power_uw);%the V parameter for each run
    pol_cont = (pol_data_table.max_power_uw-pol_data_table.min_power_uw)...
        ./(pol_data_table.max_power_uw+pol_data_table.min_power_uw);%contrast
    pol_theta = pol_data_table.min_power_angle_deg.*pi/180; %using min pow angle

    if min(pol_cont)<0
        warning('pol data has produced a contrast value less than zero')
    end
elseif strcmp(pol_opts.location,'pre_left')
    % data 2019-03-20 in lab book
    pol_data_table_pre=readtable('.\data\polz_data_pre_window_left.csv');
     %pol_data_val.Properties.VariableNames
    %'qwp_angle_deg,hwp_angle_deg,max_power_uw,max_power_angle_deg,min_power_uw,min_power_angle_deg,handedness,'
    pol_data_table_pre.min_power_angle_deg=pre_angle_flip*(pol_data_table_pre.min_power_angle_deg+pre_rotation_angle);
    pol_data_table_pre.handedness=pol_data_table_pre.handedness*pre_hand_flip;
    mean_power=mean(pol_data_table_pre.max_power_uw+pol_data_table_pre.min_power_uw);
    fprintf('mean power in pre_left %f\n',mean_power)
    
    pol_data_table_post=readtable('.\data\polz_data_post_window.csv');
    pol_data_table_post.min_power_angle_deg=pol_data_table_post.min_power_angle_deg+post_rotation_angle;
    %to keep things in a consistent reference frame rotate the angles so that they are relative to up
    %todo, do for all the other angles
    % use the qwp data from this
    mask_post_data_to_use=~isnan(pol_data_table_post.qwp_angle_deg);
    pol_data_table_post=pol_data_table_post(mask_post_data_to_use,:);
    
    commom_var_names=intersect(pol_data_table_post.Properties.VariableNames,...
                                pol_data_table_pre.Properties.VariableNames);
    
    pol_data_table=[pol_data_table_pre(:,commom_var_names);pol_data_table_post(:,commom_var_names)];
    
    
    %if we want to use the observation method
    pol_v =pol_data_table.handedness.*2.*sqrt(pol_data_table.max_power_uw.*pol_data_table.min_power_uw)...
        ./(pol_data_table.max_power_uw+pol_data_table.min_power_uw);%the V parameter for each run
    pol_cont = (pol_data_table.max_power_uw-pol_data_table.min_power_uw)...
        ./(pol_data_table.max_power_uw+pol_data_table.min_power_uw);%contrast
    pol_theta = pol_data_table.min_power_angle_deg.*pi/180; %using min pow angle

    if min(pol_cont)<0
        warning('pol data has produced a contrast value less than zero')
    end
elseif strcmp(pol_opts.location,'pre_right')
    % data 2019-03-20 in lab book
    pol_data_table_pre=readtable('.\data\polz_data_pre_window_right.csv');
    %pol_data_val.Properties.VariableNames
    %'qwp_angle_deg,hwp_angle_deg,max_power_uw,max_power_angle_deg,min_power_uw,min_power_angle_deg,handedness,'
    pol_data_table_pre.min_power_angle_deg=pre_angle_flip*(pol_data_table_pre.min_power_angle_deg+pre_rotation_angle);
    pol_data_table_pre.handedness=pol_data_table_pre.handedness*pre_hand_flip;
    mean_power=mean(pol_data_table_pre.max_power_uw+pol_data_table_pre.min_power_uw);
    fprintf('mean power in pre_right %f\n',mean_power)
    
    pol_data_table_post=readtable('.\data\polz_data_post_window.csv');
    pol_data_table_post.min_power_angle_deg=pol_data_table_post.min_power_angle_deg+post_rotation_angle;
    %to keep things in a consistent reference frame rotate the angles so that they are relative to up
    %todo, do for all the other angles
    % use the qwp data from this
    mask_post_data_to_use=~isnan(pol_data_table_post.qwp_angle_deg);
    pol_data_table_post=pol_data_table_post(mask_post_data_to_use,:);
    
    commom_var_names=intersect(pol_data_table_post.Properties.VariableNames,...
                                pol_data_table_pre.Properties.VariableNames);
    
    pol_data_table=[pol_data_table_pre(:,commom_var_names);pol_data_table_post(:,commom_var_names)];
    
    
    %if we want to use the observation method
    pol_v =pol_data_table.handedness.*2.*sqrt(pol_data_table.max_power_uw.*pol_data_table.min_power_uw)...
        ./(pol_data_table.max_power_uw+pol_data_table.min_power_uw);%the V parameter for each run
    pol_cont = (pol_data_table.max_power_uw-pol_data_table.min_power_uw)...
        ./(pol_data_table.max_power_uw+pol_data_table.min_power_uw);%contrast
    pol_theta = pol_data_table.min_power_angle_deg.*pi/180; %using min pow angle

    if min(pol_cont)<0
        warning('pol data has produced a contrast value less than zero')
    end
end
    

if strcmp(pol_opts.predict_method,'only_data')
    
    iimax=numel(pol_opts.hwp);
    for ii=1:iimax
        qwp_val_query=pol_opts.qwp(ii);
        %get mask treating nans as equal
        is_qwp_eq= (pol_data_table.qwp_angle_deg==qwp_val_query) | (isnan(qwp_val_query) & isnan(pol_data_table.qwp_angle_deg));
        hwp_val_query=pol_opts.hwp(ii);
        is_hwp_eq= (pol_data_table.hwp_angle_deg==hwp_val_query) | (isnan(hwp_val_query) & isnan(pol_data_table.hwp_angle_deg));
        query_match_data= is_qwp_eq & is_hwp_eq;
        if sum(query_match_data)>0
            data_match_idx=find(query_match_data);
            if numel(data_match_idx)>1
                warning('multiple data matches found for query point, using first match')
                data_match_idx=data_match_idx(1);
            end
            %warning('using direct data for unmodeled points, may be missing modulo')
            out_polz_state.cont.val(ii)=bound(pol_cont(data_match_idx),0,1);
            out_polz_state.theta.val(ii)=pol_theta(data_match_idx);%mod(,pi);
            out_polz_state.v.val(ii)=pol_v(data_match_idx);
        end
    end

elseif sum(strcmp(pol_opts.predict_method,{'full_fit_pref_fit','full_fit_pref_data','full_fit_only'}))>0
    % make a model for the case of having qwp+hwp, and just qwp
    
    mask_qwp_and_hwp=~isnan(pol_data_table.qwp_angle_deg) & ~isnan(pol_data_table.hwp_angle_deg);
    hwp_const_val= mode(pol_data_table.hwp_angle_deg(~isnan(pol_data_table.qwp_angle_deg)));
    mask_qwp_and_hwp_const=~isnan(pol_data_table.qwp_angle_deg) & pol_data_table.hwp_angle_deg==hwp_const_val;
    mask_hwp_only=~isnan(pol_data_table.hwp_angle_deg) & isnan(pol_data_table.qwp_angle_deg);
    mask_qwp_only=isnan(pol_data_table.hwp_angle_deg) & ~isnan(pol_data_table.qwp_angle_deg);
    stfig('polz fits');
    clf
    if sum(mask_qwp_only)>0
        warning('not set up to deal with just qwp')
    end
    
    %% fit the 4th stokes parameter for the hwp data
    beta0 = [-0.18297,1.5186,-0.033813];
    fit_V_hwp_only = fitnlm(pol_data_table.hwp_angle_deg(mask_hwp_only),pol_v(mask_hwp_only),sin_mdl_fixed_freq4,beta0);
    subplot(2,3,1)
    plot(pol_data_table.hwp_angle_deg(mask_hwp_only),pol_v(mask_hwp_only),'x')
    hold on
    hwp_samp_vec=col_vec(linspace(min(pol_data_table.hwp_angle_deg(mask_hwp_only)),max(pol_data_table.hwp_angle_deg(mask_hwp_only)),1e3));
    [v_samp_fit_val,V_samp_fit_ci]=predict(fit_V_hwp_only,hwp_samp_vec);
    plot(hwp_samp_vec,v_samp_fit_val,'-k')
    plot(hwp_samp_vec,V_samp_fit_ci,'-b')
    plot(hwp_samp_vec,sin_mdl_fixed_freq4(beta0,hwp_samp_vec),'-r')
    hold off
    xlabel('hwp angle')
    ylabel('4th stokes')
    title('hwp only')
        
    %fit_theta_hwp = fitnlm(mod(hwp_ang+21,90),mod(theta(1:37),pi),lin_mdl,beta0);
    %beta0 = [1.0,2.0139,-0.41317,-0.023839];
    %fit_d_p_qwp = fitnlm(qwp_ang,d_p(38:63),sin_mdl_abs,beta0);
    
    %% fit the theta, 2nd,3rd stokes parameter angle for the hwp data
    hwp_wraped=pol_data_table.hwp_angle_deg(mask_hwp_only);
    hwp_wraped=mod(hwp_wraped,90);
    theta_wraped=pol_theta(mask_hwp_only);
    theta_wraped=mod(theta_wraped,pi);
    [~,sort_idx]=sort(hwp_wraped);
    theta_wraped(sort_idx)=unwrap(theta_wraped(sort_idx)*2)/2;
    
    beta0=[0.034,0.5];
    fit_theta_hwp_only = fitnlm(hwp_wraped,theta_wraped,lin_mdl,beta0);
    subplot(2,3,2)
    plot(hwp_wraped,theta_wraped,'x')
    hold on
    hwp_samp_vec=col_vec(linspace(min(hwp_wraped),max(hwp_wraped),1e3));
    [theta_samp_fit_val,theta_samp_fit_ci]=predict(fit_theta_hwp_only,hwp_samp_vec);
    plot(hwp_samp_vec,theta_samp_fit_val,'-k')
    plot(hwp_samp_vec,theta_samp_fit_ci,'-b')
    plot(hwp_samp_vec,lin_mdl(beta0,hwp_samp_vec),'-r')
    hold off
    xlabel('hwp angle')
    ylabel('theta')
    title('hwp only')
    
    %% fit the theta 2nd,3rd stokes parameter contrast for the hwp data

    poly_mdl=@(b,x) taylor_series(x,b,90/2);
    %beta0 = [1,0,-0.001,0.001,1e-6,1e-6,1e-7,1e-8];
    beta0=[9.923421e-01,3.107072e-04,5.344973e-05,-8.093548e-06,-1.188569e-06,1.759216e-07,1.218760e-08,-2.301748e-09];
     opt = statset('TolFun',1e-10,'TolX',1e-10,...
            'MaxIter',1e4,... %1e4
            'UseParallel',1);
    %coef_names={'amp','phase','offset'};
    fit_cont_hwp_only = fitnlm(hwp_wraped,pol_cont(mask_hwp_only),poly_mdl,beta0,...
        'options',opt);         %'CoefficientNames',coef_names,...
    subplot(2,3,3)
    plot(hwp_wraped,pol_cont(mask_hwp_only),'x')
    hold on
    hwp_samp_vec=col_vec(linspace(min(hwp_wraped),max(hwp_wraped),1e3));
    [cont_samp_fit_val,cont_samp_fit_ci]=predict(fit_cont_hwp_only,hwp_samp_vec);
%     theta_samp_fit_val=bound(theta_samp_fit_val,0,1);
    plot(hwp_samp_vec,cont_samp_fit_val,'-k')
    plot(hwp_samp_vec,cont_samp_fit_ci,'-b')
    %plot(hwp_samp_vec,poly_mdl(beta0,hwp_samp_vec),'-r')
    hold off
    xlabel('hwp angle')
    ylabel('contrast')
    title('hwp only')
    
    %% fit the 4th stokes parameter for the qwp data (with hwp=333)
     
    % fit only the dominant hwp angle here
    %idx_qwp_and_hwp_333
    beta0 = [-1.0002,-0.41317,-0.023839];
    fit_V_qwp_and_hwp_const = fitnlm(pol_data_table.qwp_angle_deg(mask_qwp_and_hwp_const),pol_v(mask_qwp_and_hwp_const),sin_mdl_fixed_freq2,beta0);
    subplot(2,3,4)
    plot(pol_data_table.qwp_angle_deg(mask_qwp_and_hwp_const),pol_v(mask_qwp_and_hwp_const),'x')
    hold on
    qwp_angle_samp_fit_hwp_const=col_vec(linspace(min(pol_data_table.qwp_angle_deg(mask_qwp_and_hwp_const)),max(pol_data_table.qwp_angle_deg(mask_qwp_and_hwp_const)),1e3));
    [v_samp_fit_val,V_samp_fit_ci]=predict(fit_V_qwp_and_hwp_const,qwp_angle_samp_fit_hwp_const);
    plot(qwp_angle_samp_fit_hwp_const,v_samp_fit_val,'-k')
    plot(qwp_angle_samp_fit_hwp_const,V_samp_fit_ci,'-b')
    plot(qwp_angle_samp_fit_hwp_const,sin_mdl_fixed_freq2(beta0,qwp_angle_samp_fit_hwp_const),'-r')
    hold off
    xlabel('qwp angle')
    ylabel('4th stokes')
    title(sprintf('qwp & hwp=%.1f',hwp_const_val))
    
    %% fit the theta 2nd,3rd stokes angle parameter for the qwp data (with hwp=333)
    qwp_wraped=mod(pol_data_table.qwp_angle_deg(mask_qwp_and_hwp_const),90);
    theta_wraped=mod(pol_theta(mask_qwp_and_hwp_const),pi);
    [~,sort_idx]=sort(qwp_wraped);
    theta_wraped(sort_idx)=unwrap(theta_wraped(sort_idx)*2)/2;
    
    beta0=[0.5,-1,2.5];
    opt = statset('TolFun',1e-10,'TolX',1e-10,...
            'MaxIter',1e4,... %1e4
            'UseParallel',1);
    %sin_mdl = @(b,x) b(1).*sin(b(2).*x.*pi/180+b(3))+b(4);
    coef_names={'amp','phase','offset'};
    fit_theta_qwp_and_hwp_const = fitnlm(pol_data_table.qwp_angle_deg(mask_qwp_and_hwp_const),theta_wraped,sin_mdl_fixed_freq4,beta0,...
        'CoefficientNames',coef_names,...
        'options',opt);
    subplot(2,3,5)
    plot(pol_data_table.qwp_angle_deg(mask_qwp_and_hwp_const),theta_wraped,'x')
    hold on
    v_hwp_angle_samp_fit_hwp_const=col_vec(linspace(min(pol_data_table.qwp_angle_deg(mask_qwp_and_hwp_const)),max(pol_data_table.qwp_angle_deg(mask_qwp_and_hwp_const)),1e3));
    [theta_samp_fit_val,theta_samp_fit_ci]=predict(fit_theta_qwp_and_hwp_const,v_hwp_angle_samp_fit_hwp_const);
    plot(v_hwp_angle_samp_fit_hwp_const,theta_samp_fit_val,'-k')
    plot(v_hwp_angle_samp_fit_hwp_const,theta_samp_fit_ci,'-b')
    plot(v_hwp_angle_samp_fit_hwp_const,sin_mdl_fixed_freq4(beta0,v_hwp_angle_samp_fit_hwp_const),'-r')
    hold off
    xlabel('qwp angle')
    ylabel('theta')
    title(sprintf('qwp & hwp=%.1f',hwp_const_val))
    
    %% fit the theta 2nd,3rd stokes parameter contrast for the qwp data (with hwp=333)
    
    beta0 = [0.3,1.5186,0.5];
    opt = statset('TolFun',1e-10,'TolX',1e-10,...
            'MaxIter',1e4,... %1e4
            'UseParallel',1);
    coef_names={'amp','phase','offset'};
    fit_cont_qwp_const_hwp = fitnlm(pol_data_table.qwp_angle_deg(mask_qwp_and_hwp_const),pol_cont(mask_qwp_and_hwp_const),sin_mdl_fixed_freq4,beta0,...
        'CoefficientNames',coef_names,...
        'options',opt);
    subplot(2,3,6)
    plot(pol_data_table.qwp_angle_deg(mask_qwp_and_hwp_const),pol_cont(mask_qwp_and_hwp_const),'x')
    hold on
    hwp_samp_vec=col_vec(linspace(min(pol_data_table.qwp_angle_deg(mask_qwp_and_hwp_const)),max(pol_data_table.qwp_angle_deg(mask_qwp_and_hwp_const)),1e3));
    [cont_samp_fit_val,cont_samp_fit_ci]=predict(fit_cont_qwp_const_hwp,hwp_samp_vec);
    % clip the min max value to be 1
%     cont_samp_fit_val=bound(cont_samp_fit_val,0,1);
    plot(hwp_samp_vec,cont_samp_fit_val,'-k')
    plot(hwp_samp_vec,cont_samp_fit_ci,'-b')
    plot(hwp_samp_vec,sin_mdl_fixed_freq4(beta0,hwp_samp_vec),'-r')
    hold off
    xlabel('qwp angle')
    ylabel('contrast')
    title(sprintf('qwp & hwp=%.1f',hwp_const_val))
    
    %% begin query fufillment  
    fulfilled_query=pol_opts.hwp*false; % has the query for this hwp,qwp combo been completed
    
    %%

    if strcmp(pol_opts.predict_method,'full_fit_pref_data')
        iimax=numel(pol_opts.hwp);
        for ii=1:iimax
            qwp_val_query=pol_opts.qwp(ii);
            %get mask treating nans as equal
            is_qwp_eq= (pol_data_table.qwp_angle_deg==qwp_val_query) | (isnan(qwp_val_query) & isnan(pol_data_table.qwp_angle_deg));
            hwp_val_query=pol_opts.hwp(ii);
            is_hwp_eq= (pol_data_table.hwp_angle_deg==hwp_val_query) | (isnan(hwp_val_query) & isnan(pol_data_table.hwp_angle_deg));
            query_match_data= is_qwp_eq & is_hwp_eq;
            if sum(query_match_data)>0
                warning('using direct data for unmodeled points, may be missing modulo')
                data_match_idx=find(query_match_data);
                if numel(data_match_idx)>1
                    warning('multiple data matches found for query point')
                    data_match_idx=data_match_idx(1);
                end
                
                out_polz_state.cont.val(ii)=bound(pol_cont(data_match_idx),0,1); %bound the contrast to be 1
                out_polz_state.theta.val(ii)=pol_theta(data_match_idx);  %%mod(,pi/2);
                out_polz_state.v.val(ii)=pol_v(data_match_idx);
                out_polz_state.v.unc(ii)=nan;
                out_polz_state.theta.unc(ii)=nan;
                out_polz_state.cont.unc(ii)=nan;
                fulfilled_query(ii)=true;
            end
        end
        
    end 
    
    
    %% now calulate at the remaining query points
    only_qwp=isnan(pol_opts.hwp) & ~isnan( pol_opts.qwp) & ~fulfilled_query;
    if sum( only_qwp)>0
        warning('not set up to query quater wp by itself')
    end
    
    %% query the model at the hwp only data pts that remain
    mask_query_no_qwp=~isnan(pol_opts.hwp) & isnan(pol_opts.qwp) & ~fulfilled_query;
    fulfilled_query=fulfilled_query & mask_query_no_qwp;
    query_hwp_wraped=pol_opts.hwp(mask_query_no_qwp);
    query_hwp_wraped=mod(query_hwp_wraped,90);
    
    [a,b]=predict(fit_V_hwp_only,pol_opts.hwp(mask_query_no_qwp),'Alpha',1-erf(1/sqrt(2)));
    out_polz_state.v.val(mask_query_no_qwp)=a;
    out_polz_state.v.unc(mask_query_no_qwp)=range(b,2)/2;
    [a,b]=predict(fit_theta_hwp_only,query_hwp_wraped,'Alpha',1-erf(1/sqrt(2)));

    out_polz_state.theta.val(mask_query_no_qwp)=mod(a,pi);
    out_polz_state.theta.unc(mask_query_no_qwp)=range(b,2)/2;
    [a,b]=predict(fit_cont_hwp_only,query_hwp_wraped,'Alpha',1-erf(1/sqrt(2)));
    out_polz_state.cont.val(mask_query_no_qwp)=a;%bound(a,0,1);
    out_polz_state.cont.unc(mask_query_no_qwp)=range(b,2)/2;
    
    subplot(2,3,1)
    hold on
    plot(pol_opts.hwp(mask_query_no_qwp), out_polz_state.v.val(mask_query_no_qwp),'ko')
    hold off
    subplot(2,3,2)
    hold on
    plot(query_hwp_wraped,out_polz_state.theta.val(mask_query_no_qwp),'ko')
    hold off
    subplot(2,3,3)
    hold on
    plot(query_hwp_wraped,out_polz_state.cont.val(mask_query_no_qwp),'ko')
    hold off
    %%

    mask_query_both_with_hwp_const=pol_opts.hwp==hwp_const_val & ~isnan(pol_opts.qwp) & ~fulfilled_query;
    fulfilled_query=fulfilled_query & mask_query_both_with_hwp_const;
    [a,b]=predict(fit_V_qwp_and_hwp_const,pol_opts.qwp(mask_query_both_with_hwp_const),'Alpha',1-erf(1/sqrt(2)));
    out_polz_state.v.val(mask_query_both_with_hwp_const)=a;
    out_polz_state.v.unc(mask_query_both_with_hwp_const)=range(b,2)/2;
    [a,b]=predict(fit_theta_qwp_and_hwp_const,pol_opts.qwp(mask_query_both_with_hwp_const),'Alpha',1-erf(1/sqrt(2)));
    out_polz_state.theta.val(mask_query_both_with_hwp_const)=mod(a,pi);
    out_polz_state.theta.unc(mask_query_both_with_hwp_const)=range(b,2)/2;
    [a,b]=predict(fit_cont_qwp_const_hwp,pol_opts.qwp(mask_query_both_with_hwp_const),'Alpha',1-erf(1/sqrt(2)));
    out_polz_state.cont.val(mask_query_both_with_hwp_const)=a;%bound(a,0,1);
    out_polz_state.cont.unc(mask_query_both_with_hwp_const)=range(b,2)/2;
    
    subplot(2,3,4)
    hold on
    plot(pol_opts.qwp(mask_query_both_with_hwp_const),  out_polz_state.v.val(mask_query_both_with_hwp_const),'ko')
    hold off
    subplot(2,3,5)
    hold on
    plot(pol_opts.qwp(mask_query_both_with_hwp_const), out_polz_state.theta.val(mask_query_both_with_hwp_const),'ko')
    hold off
    subplot(2,3,6)
    hold on
    plot(pol_opts.qwp(mask_query_both_with_hwp_const),out_polz_state.cont.val(mask_query_both_with_hwp_const),'ko')
    hold off

    %% try to see if the rest of the query points can just be extracted from the data
    if  strcmp(pol_opts.predict_method,'full_fit_pref_fit')
        %mask_both_with_hwp_notconst=pol_opts.hwp~=hwp_const_val & ~isnan(pol_opts.hwp) & ~isnan(pol_opts.qwp);
        mask_unmatched=~mask_query_both_with_hwp_const & ~mask_query_no_qwp & ~fulfilled_query;
        idxs_unmatched=find(mask_unmatched);
        iimax=sum(mask_unmatched);
        for ii=1:iimax
            query_idx=idxs_unmatched(ii);
            qwp_val_query=pol_opts.qwp(query_idx);
            %get mask treating nans as equal
            is_qwp_eq= (pol_data_table.qwp_angle_deg==qwp_val_query) | (isnan(qwp_val_query) & isnan(pol_data_table.qwp_angle_deg));
            hwp_val_query=pol_opts.hwp(query_idx);
            is_hwp_eq= (pol_data_table.hwp_angle_deg==hwp_val_query) | (isnan(hwp_val_query) & isnan(pol_data_table.hwp_angle_deg));
            query_match_data= is_qwp_eq & is_hwp_eq;
            if sum(query_match_data)>0
                warning('using direct data for unmodeled points, may be missing modulo')
                data_match_idx=find(query_match_data);
                if numel(data_match_idx)>1
                    warning('multiple data matches found for query point')
                    data_match_idx=data_match_idx(1);
                end
                
                out_polz_state.cont.val(query_idx)=bound(pol_cont(data_match_idx),0,1);
                out_polz_state.theta.val(query_idx)=mod(pol_theta(data_match_idx),pi);
                out_polz_state.v.val(query_idx)=pol_v(data_match_idx);

                out_polz_state.v.unc(query_idx)=nan;
                out_polz_state.theta.unc(query_idx)=nan;
                out_polz_state.cont.unc(query_idx)=nan;
            end
        end
    end
    
    
 
    
    
    
    
    
    
elseif sum(strcmp(pol_opts.predict_method,{'interp_only','interp_pref_data','interp_pref_interp'}))>0
     % make a model for the case of having qwp+hwp, and just qwp
    
    mask_qwp_and_hwp=~isnan(pol_data_table.qwp_angle_deg) & ~isnan(pol_data_table.hwp_angle_deg);
    hwp_const_val= mode(pol_data_table.hwp_angle_deg(~isnan(pol_data_table.qwp_angle_deg)));
    mask_qwp_and_hwp_const=~isnan(pol_data_table.qwp_angle_deg) & pol_data_table.hwp_angle_deg==hwp_const_val;
    mask_hwp_only=~isnan(pol_data_table.hwp_angle_deg) & isnan(pol_data_table.qwp_angle_deg);
    mask_qwp_only=isnan(pol_data_table.hwp_angle_deg) & ~isnan(pol_data_table.qwp_angle_deg);
    stfig('polz fits');
    clf
    if sum(mask_qwp_only)>0
        warning('not set up to deal with just qwp')
    end
    
    %% fit the 4th stokes parameter for the hwp data
    hwp_raw=pol_data_table.hwp_angle_deg(mask_hwp_only);
    if pol_opts.wrap_hwp
        hwp_wraped=mod(hwp_raw,90);
    else
        hwp_wraped=hwp_raw;
    end
    
    
    subplot(2,3,1)
    plot(hwp_wraped,pol_v(mask_hwp_only),'x')
    hold on
    %annoyingly we need to consolidate points in x,y
    [x_hwp_v,y_hwp_v]=consolidator11(hwp_wraped,pol_v(mask_hwp_only));
    hwp_samp_vec=col_vec(linspace(min(hwp_wraped),max(hwp_wraped),1e3));
    v_samp_fit_val=pchip(x_hwp_v,y_hwp_v,hwp_samp_vec);
    plot(hwp_samp_vec,v_samp_fit_val,'-k')
    hold off
    xlabel('hwp angle')
    ylabel('4th stokes')
    title('hwp only')
        
    %fit_theta_hwp = fitnlm(mod(hwp_ang+21,90),mod(theta(1:37),pi),lin_mdl,beta0);
    %beta0 = [1.0,2.0139,-0.41317,-0.023839];
    %fit_d_p_qwp = fitnlm(qwp_ang,d_p(38:63),sin_mdl_abs,beta0);
    
    %% fit the theta, 2nd,3rd stokes parameter angle for the hwp data
    theta_wraped_hwp=pol_theta(mask_hwp_only);
    theta_wraped_hwp=mod(theta_wraped_hwp,pi);
    [~,sort_idx]=sort(hwp_wraped);
    theta_wraped_hwp(sort_idx)=unwrap(theta_wraped_hwp(sort_idx)*2)/2;
    
    subplot(2,3,2)
    plot(hwp_wraped,theta_wraped_hwp,'x')
    hold on
    [x_hwp_theta,y_hwp_theta]=consolidator11(hwp_wraped,theta_wraped_hwp);
    hwp_samp_vec=col_vec(linspace(min(hwp_wraped),max(hwp_wraped),1e3));
    theta_samp_fit_val=pchip(x_hwp_theta,y_hwp_theta,hwp_samp_vec);
    plot(hwp_samp_vec,theta_samp_fit_val,'-k')
    hold off
    xlabel('hwp angle')
    ylabel('theta')
    title('hwp only')
    
    %% fit the theta 2nd,3rd stokes parameter contrast for the hwp data

    subplot(2,3,3)
    plot(hwp_wraped,pol_cont(mask_hwp_only),'x')
    hold on
    [x_hwp_cont,y_hwp_cont]=consolidator11(hwp_wraped,pol_cont(mask_hwp_only));
    hwp_samp_vec=col_vec(linspace(min(hwp_wraped),max(hwp_wraped),1e3));
    theta_samp_fit_val=pchip(x_hwp_cont,y_hwp_cont,hwp_samp_vec);
%     theta_samp_fit_val=bound(theta_samp_fit_val,0,1);
    plot(hwp_samp_vec,theta_samp_fit_val,'-k')
    hold off
    xlabel('hwp angle')
    ylabel('contrast')
    title('hwp only')
    
    %% fit the 4th stokes parameter for the qwp data (with hwp=333)
     
    % fit only the dominant hwp angle here
    %idx_qwp_and_hwp_333

    subplot(2,3,4)
    plot(pol_data_table.qwp_angle_deg(mask_qwp_and_hwp_const),pol_v(mask_qwp_and_hwp_const),'x')
    hold on
    [x_qwp_v,y_qwp_v]=consolidator11(pol_data_table.qwp_angle_deg(mask_qwp_and_hwp_const),pol_v(mask_qwp_and_hwp_const));
    qwp_angle_samp_fit_hwp_const=col_vec(linspace(min(pol_data_table.qwp_angle_deg(mask_qwp_and_hwp_const)),max(pol_data_table.qwp_angle_deg(mask_qwp_and_hwp_const)),1e3));
    v_samp_fit_val=pchip(x_qwp_v,y_qwp_v,qwp_angle_samp_fit_hwp_const);
    plot(qwp_angle_samp_fit_hwp_const,v_samp_fit_val,'-k')
    hold off
    xlabel('qwp angle')
    ylabel('4th stokes')
    title(sprintf('qwp & hwp=%.1f',hwp_const_val))
    
    %% fit the theta 2nd,3rd stokes angle parameter for the qwp data (with hwp=333)
    qwp_wraped=mod(pol_data_table.qwp_angle_deg(mask_qwp_and_hwp_const),90);
    theta_wraped=mod(pol_theta(mask_qwp_and_hwp_const),pi);
    [~,sort_idx]=sort(qwp_wraped);
    theta_wraped(sort_idx)=unwrap(theta_wraped(sort_idx)*2)/2;
    subplot(2,3,5)
    plot(pol_data_table.qwp_angle_deg(mask_qwp_and_hwp_const),theta_wraped,'x')
    hold on
    [x_qwp_theta,y_qwp_theta]=consolidator11(pol_data_table.qwp_angle_deg(mask_qwp_and_hwp_const),theta_wraped);
    hwp_angle_samp_fit_hwp_const=col_vec(linspace(min(pol_data_table.qwp_angle_deg(mask_qwp_and_hwp_const)),max(pol_data_table.qwp_angle_deg(mask_qwp_and_hwp_const)),1e3));
    theta_samp_fit_val=pchip(x_qwp_theta,y_qwp_theta,hwp_angle_samp_fit_hwp_const);
    plot(hwp_angle_samp_fit_hwp_const,theta_samp_fit_val,'-k')
    hold off
    xlabel('qwp angle')
    ylabel('theta')
    title(sprintf('qwp & hwp=%.1f',hwp_const_val))
    
    %% fit the contrast 2nd,3rd stokes parameter contrast for the qwp data (with hwp=333)
    
    subplot(2,3,6)
    plot(pol_data_table.qwp_angle_deg(mask_qwp_and_hwp_const),pol_cont(mask_qwp_and_hwp_const),'x')
    hold on
    [x_qwp_cont,y_qwp_cont]=consolidator11(pol_data_table.qwp_angle_deg(mask_qwp_and_hwp_const),pol_cont(mask_qwp_and_hwp_const));
    contrast_samp_fit_val=pchip(x_qwp_cont,y_qwp_cont,hwp_angle_samp_fit_hwp_const);
    % clip the min max value to be 1
%     contrast_samp_fit_val=bound(contrast_samp_fit_val,0,1);
    plot(hwp_angle_samp_fit_hwp_const,contrast_samp_fit_val,'-k')
    hold off
    xlabel('qwp angle')
    ylabel('contrast')
    title(sprintf('qwp & hwp=%.1f',hwp_const_val))
    
    %% begin query fufillment  
    fulfilled_query=pol_opts.hwp*false; % has the query for this hwp,qwp combo been completed
    

    if strcmp(pol_opts.predict_method,'interp_pref_data')
        iimax=numel(pol_opts.hwp);
        for ii=1:iimax
            qwp_val_query=pol_opts.qwp(ii);
            %get mask treating nans as equal
            is_qwp_eq= (pol_data_table.qwp_angle_deg==qwp_val_query) | (isnan(qwp_val_query) & isnan(pol_data_table.qwp_angle_deg));
            hwp_val_query=pol_opts.hwp(ii);
            is_hwp_eq= (pol_data_table.hwp_angle_deg==hwp_val_query) | (isnan(hwp_val_query) & isnan(pol_data_table.hwp_angle_deg));
            query_match_data= is_qwp_eq & is_hwp_eq;
            if sum(query_match_data)>0
                warning('using direct data for unmodeled points, may be missing modulo')
                data_match_idx=find(query_match_data);
                if numel(data_match_idx)>1
                    warning('multiple data matches found for query point')
                    data_match_idx=data_match_idx(1);
                end
                
                out_polz_state.cont.val(ii)=bound(pol_cont(data_match_idx),0,1); %bound the contrast to be 1
                out_polz_state.theta.val(ii)=mod(pol_theta(data_match_idx),2);
                out_polz_state.v.val(ii)=pol_v(data_match_idx);
                out_polz_state.v.unc(ii)=nan;
                out_polz_state.theta.unc(ii)=nan;
                out_polz_state.cont.unc(ii)=nan;
                fulfilled_query(ii)=true;
            end
        end
        
    end 
    
    
    %% now calulate at the remaining query points
    only_qwp=isnan(pol_opts.hwp) & ~isnan( pol_opts.qwp) & ~fulfilled_query;
    if sum( only_qwp)>0
        warning('not set up to query quater wp by itself')
        out_theta_v(only_qwp,1)=nan;
        out_theta_v(only_qwp,2)=nan;
    end
    
    %% query the model at the hwp only data pts that remain
    mask_query_no_qwp=~isnan(pol_opts.hwp) & isnan(pol_opts.qwp) & ~fulfilled_query;
    query_hwp_raw=pol_opts.hwp(mask_query_no_qwp);
    if pol_opts.wrap_hwp
        query_hwp_wraped=mod(query_hwp_raw,90);
    else
        query_hwp_wraped=query_hwp_raw;
    end
    
    out_polz_state.v.val(mask_query_no_qwp)=pchip(x_hwp_v,y_hwp_v,query_hwp_wraped);
    out_polz_state.theta.val(mask_query_no_qwp)=pchip(x_hwp_theta,y_hwp_theta,query_hwp_wraped);
    a=pchip(x_hwp_cont,y_hwp_cont,query_hwp_wraped);
    out_polz_state.cont.val(mask_query_no_qwp)=a;%bound(a,0,1);
    
    
    subplot(2,3,1)
    hold on
    plot(query_hwp_wraped, out_polz_state.v.val(mask_query_no_qwp),'ko')
    hold off
    subplot(2,3,2)
    hold on
    plot(query_hwp_wraped,out_polz_state.theta.val(mask_query_no_qwp),'ko')
    hold off
    subplot(2,3,3)
    hold on
    plot(query_hwp_wraped,out_polz_state.cont.val(mask_query_no_qwp),'ko')
    hold off
    %%

    mask_query_both_with_hwp_const=pol_opts.hwp==hwp_const_val & ~isnan(pol_opts.qwp) & ~fulfilled_query;

    out_polz_state.v.val(mask_query_both_with_hwp_const)=pchip(x_qwp_v,y_qwp_v,pol_opts.qwp(mask_query_both_with_hwp_const));
    out_polz_state.theta.val(mask_query_both_with_hwp_const)=pchip(x_qwp_theta,y_qwp_theta,pol_opts.qwp(mask_query_both_with_hwp_const));
    a=pchip(x_qwp_cont,y_qwp_cont,pol_opts.qwp(mask_query_both_with_hwp_const));
    out_polz_state.cont.val(mask_query_both_with_hwp_const)=a;%bound(a,0,1);
    
    subplot(2,3,4)
    hold on
    plot(pol_opts.qwp(mask_query_both_with_hwp_const),  out_polz_state.v.val(mask_query_both_with_hwp_const),'ko')
    hold off
    subplot(2,3,5)
    hold on
    plot(pol_opts.qwp(mask_query_both_with_hwp_const), out_polz_state.theta.val(mask_query_both_with_hwp_const),'ko')
    hold off
    subplot(2,3,6)
    hold on
    plot(pol_opts.qwp(mask_query_both_with_hwp_const),out_polz_state.cont.val(mask_query_both_with_hwp_const),'ko')
    hold off

    %% try to see if the rest of the query points can just be extracted from the data
    if  strcmp(pol_opts.predict_method,'interp_pref_interp')
        %mask_both_with_hwp_notconst=pol_opts.hwp~=hwp_const_val & ~isnan(pol_opts.hwp) & ~isnan(pol_opts.qwp);
        mask_unmatched=~mask_query_both_with_hwp_const & ~mask_query_no_qwp & ~fulfilled_query;
        idxs_unmatched=find(mask_unmatched);
        iimax=sum(mask_unmatched);
        for ii=1:iimax
            query_idx=idxs_unmatched(ii);
            qwp_val_query=pol_opts.qwp(query_idx);
            %get mask treating nans as equal
            is_qwp_eq= (pol_data_table.qwp_angle_deg==qwp_val_query) | (isnan(qwp_val_query) & isnan(pol_data_table.qwp_angle_deg));
            hwp_val_query=pol_opts.hwp(query_idx);
            is_hwp_eq= (pol_data_table.hwp_angle_deg==hwp_val_query) | (isnan(hwp_val_query) & isnan(pol_data_table.hwp_angle_deg));
            query_match_data= is_qwp_eq & is_hwp_eq;
            if sum(query_match_data)>0
                warning('using direct data for unmodeled points, may be missing modulo')
                data_match_idx=find(query_match_data);
                if numel(data_match_idx)>1
                    warning('multiple data matches found for query point')
                    data_match_idx=data_match_idx(1);
                end
                
                out_polz_state.cont.val(query_idx)=bound(pol_cont(data_match_idx),0,1);
                out_polz_state.theta.val(query_idx)=mod(pol_theta(data_match_idx),pi);
                out_polz_state.v.val(query_idx)=pol_v(data_match_idx);

                out_polz_state.v.unc(query_idx)=nan;
                out_polz_state.theta.unc(query_idx)=nan;
                out_polz_state.cont.unc(query_idx)=nan;
            end
        end
    end

    
elseif sum(strcmp(pol_opts.predict_method,{'gauss_only','gauss_pref_data','gauss_pref_interp'}))>0
    % make a model for the case of having qwp+hwp, and just qwp
    
    mask_qwp_and_hwp=~isnan(pol_data_table.qwp_angle_deg) & ~isnan(pol_data_table.hwp_angle_deg);
    hwp_const_val= mode(pol_data_table.hwp_angle_deg(~isnan(pol_data_table.qwp_angle_deg)));
    mask_qwp_and_hwp_const=~isnan(pol_data_table.qwp_angle_deg) & pol_data_table.hwp_angle_deg==hwp_const_val;
    mask_hwp_only=~isnan(pol_data_table.hwp_angle_deg) & isnan(pol_data_table.qwp_angle_deg);
    mask_qwp_only=isnan(pol_data_table.hwp_angle_deg) & ~isnan(pol_data_table.qwp_angle_deg);
    stfig('polz fits');
    clf
    if sum(mask_qwp_only)>0
        warning('not set up to deal with just qwp')
    end
    
    %% fit the 4th stokes parameter for the hwp data
    hwp_raw=pol_data_table.hwp_angle_deg(mask_hwp_only);
    if pol_opts.wrap_hwp
        hwp_wraped=mod(hwp_raw,90);
    else
        hwp_wraped=hwp_raw;
    end
    
    subplot(2,3,1)
    plot(hwp_wraped,pol_v(mask_hwp_only),'x')
    hold on

    
    x_hwp_v=hwp_wraped;
    y_hwp_v=pol_v(mask_hwp_only);
    hwp_samp_vec=col_vec(linspace(min(hwp_wraped),max(hwp_wraped),1e3));
    sigma_hwp_v=pol_opts.smoothing;
    v_samp_fit_val=gauss_weighted_interp(x_hwp_v,y_hwp_v,hwp_samp_vec,sigma_hwp_v);
    
    plot(hwp_samp_vec,v_samp_fit_val,'-k')
    hold off
    xlabel('hwp angle')
    ylabel('4th stokes')
    title('hwp only')
        
    %% fit the theta, 2nd,3rd stokes parameter angle for the hwp data
    theta_wraped_hwp=pol_theta(mask_hwp_only);
    theta_wraped_hwp=mod(theta_wraped_hwp,pi);
    [~,sort_idx]=sort(hwp_wraped);
    theta_wraped_hwp(sort_idx)=unwrap(theta_wraped_hwp(sort_idx)*2)/2;
    
    subplot(2,3,2)
    plot(hwp_wraped,theta_wraped_hwp,'x')
    hold on
    x_hwp_theta=hwp_wraped;
    y_hwp_theta=theta_wraped_hwp;
    hwp_samp_vec=col_vec(linspace(min(hwp_wraped),max(hwp_wraped),1e3));
    sigma_hwp_theta=pol_opts.smoothing;
    theta_samp_fit_val=gauss_weighted_interp(x_hwp_theta,y_hwp_theta,hwp_samp_vec,sigma_hwp_theta);
    plot(hwp_samp_vec,theta_samp_fit_val,'-k')
    hold off
    xlabel('hwp angle')
    ylabel('theta')
    title('hwp only')
    
    %% fit the theta 2nd,3rd stokes parameter contrast for the hwp data

    subplot(2,3,3)
    plot(hwp_wraped,pol_cont(mask_hwp_only),'x')
    hold on
    x_hwp_cont=hwp_wraped;
    y_hwp_cont=pol_cont(mask_hwp_only);
    hwp_samp_vec=col_vec(linspace(min(hwp_wraped),max(hwp_wraped),1e3));
    sigma_hwp_cont=pol_opts.smoothing;
    theta_samp_fit_val=gauss_weighted_interp(x_hwp_cont,y_hwp_cont,hwp_samp_vec,sigma_hwp_cont);
    theta_samp_fit_val=bound(theta_samp_fit_val,0,1);
    plot(hwp_samp_vec,theta_samp_fit_val,'-k')
    hold off
    xlabel('hwp angle')
    ylabel('contrast')
    title('hwp only')
    
    %% fit the 4th stokes parameter for the qwp data (with hwp=333)
     
    % fit only the dominant hwp angle here
    %idx_qwp_and_hwp_333

    subplot(2,3,4)
    plot(pol_data_table.qwp_angle_deg(mask_qwp_and_hwp_const),pol_v(mask_qwp_and_hwp_const),'x')
    hold on
    x_qwp_v=pol_data_table.qwp_angle_deg(mask_qwp_and_hwp_const);
    y_qwp_v=pol_v(mask_qwp_and_hwp_const);
    qwp_angle_samp_fit_hwp_const=col_vec(linspace(min(pol_data_table.qwp_angle_deg(mask_qwp_and_hwp_const)),max(pol_data_table.qwp_angle_deg(mask_qwp_and_hwp_const)),1e3));
    sigma_qwp_v=pol_opts.smoothing;
    v_samp_fit_val=gauss_weighted_interp(x_qwp_v,y_qwp_v,qwp_angle_samp_fit_hwp_const,sigma_qwp_v);
    plot(qwp_angle_samp_fit_hwp_const,v_samp_fit_val,'-k')
    hold off
    xlabel('qwp angle')
    ylabel('4th stokes')
    title(sprintf('qwp & hwp=%.1f',hwp_const_val))
    
    %% fit the theta 2nd,3rd stokes angle parameter for the qwp data (with hwp=333)
    qwp_wraped=mod(pol_data_table.qwp_angle_deg(mask_qwp_and_hwp_const),90);
    theta_wraped=mod(pol_theta(mask_qwp_and_hwp_const),pi);
    [~,sort_idx]=sort(qwp_wraped);
    theta_wraped(sort_idx)=unwrap(theta_wraped(sort_idx)*2)/2;
    subplot(2,3,5)
    plot(pol_data_table.qwp_angle_deg(mask_qwp_and_hwp_const),theta_wraped,'x')
    hold on
    x_qwp_theta=pol_data_table.qwp_angle_deg(mask_qwp_and_hwp_const);
    y_qwp_theta=theta_wraped;
    hwp_angle_samp_fit_hwp_const=col_vec(linspace(min(pol_data_table.qwp_angle_deg(mask_qwp_and_hwp_const)),max(pol_data_table.qwp_angle_deg(mask_qwp_and_hwp_const)),1e3));
    sigma_qwp_theta=pol_opts.smoothing;
    theta_samp_fit_val=gauss_weighted_interp(x_qwp_theta,y_qwp_theta,hwp_angle_samp_fit_hwp_const,sigma_qwp_theta);
    plot(hwp_angle_samp_fit_hwp_const,theta_samp_fit_val,'-k')
    hold off
    xlabel('qwp angle')
    ylabel('theta')
    title(sprintf('qwp & hwp=%.1f',hwp_const_val))
    
    %% fit the contrast 2nd,3rd stokes parameter contrast for the qwp data (with hwp=333)
    
    subplot(2,3,6)
    plot(pol_data_table.qwp_angle_deg(mask_qwp_and_hwp_const),pol_cont(mask_qwp_and_hwp_const),'x')
    hold on
    x_qwp_cont=pol_data_table.qwp_angle_deg(mask_qwp_and_hwp_const);
    y_qwp_cont=pol_cont(mask_qwp_and_hwp_const);
    sigma_qwp_cont=pol_opts.smoothing;
    contrast_samp_fit_val=gauss_weighted_interp(x_qwp_cont,y_qwp_cont,hwp_angle_samp_fit_hwp_const,sigma_qwp_cont);
    % clip the min max value to be 1
    contrast_samp_fit_val=bound(contrast_samp_fit_val,0,1);
    plot(hwp_angle_samp_fit_hwp_const,contrast_samp_fit_val,'-k')
    hold off
    xlabel('qwp angle')
    ylabel('contrast')
    title(sprintf('qwp & hwp=%.1f',hwp_const_val))
    
    %% begin query fufillment  
    fulfilled_query=pol_opts.hwp*false; % has the query for this hwp,qwp combo been completed
    

    if strcmp(pol_opts.predict_method,'gauss_pref_data')
        iimax=numel(pol_opts.hwp);
        for ii=1:iimax
            qwp_val_query=pol_opts.qwp(ii);
            %get mask treating nans as equal
            is_qwp_eq= (pol_data_table.qwp_angle_deg==qwp_val_query) | (isnan(qwp_val_query) & isnan(pol_data_table.qwp_angle_deg));
            hwp_val_query=pol_opts.hwp(ii);
            is_hwp_eq= (pol_data_table.hwp_angle_deg==hwp_val_query) | (isnan(hwp_val_query) & isnan(pol_data_table.hwp_angle_deg));
            query_match_data= is_qwp_eq & is_hwp_eq;
            if sum(query_match_data)>0
                warning('using direct data for unmodeled points, may be missing modulo')
                data_match_idx=find(query_match_data);
                if numel(data_match_idx)>1
                    warning('multiple data matches found for query point')
                    data_match_idx=data_match_idx(1);
                end
                
                out_polz_state.cont.val(ii)=bound(pol_cont(data_match_idx),0,1); %bound the contrast to be 1
                out_polz_state.theta.val(ii)=mod(pol_theta(data_match_idx),2);
                out_polz_state.v.val(ii)=pol_v(data_match_idx);
                out_polz_state.v.unc(ii)=nan;
                out_polz_state.theta.unc(ii)=nan;
                out_polz_state.cont.unc(ii)=nan;
                fulfilled_query(ii)=true;
            end
        end
        
    end 
    
    
    %% now calulate at the remaining query points
    only_qwp=isnan(pol_opts.hwp) & ~isnan( pol_opts.qwp) & ~fulfilled_query;
    if sum( only_qwp)>0
        warning('not set up to query quater wp by itself')
        out_theta_v(only_qwp,1)=nan;
        out_theta_v(only_qwp,2)=nan;
    end
    
    %% query the model at the hwp only data pts that remain
    mask_query_no_qwp=~isnan(pol_opts.hwp) & isnan(pol_opts.qwp) & ~fulfilled_query;
    query_hwp_raw=pol_opts.hwp(mask_query_no_qwp);
    if pol_opts.wrap_hwp
        query_hwp_wraped=mod(query_hwp_raw,90);
    else
        query_hwp_wraped=query_hwp_raw;
    end
    
    out_polz_state.v.val(mask_query_no_qwp)=gauss_weighted_interp(x_hwp_v,y_hwp_v,query_hwp_wraped,sigma_hwp_v);
    out_polz_state.theta.val(mask_query_no_qwp)=gauss_weighted_interp(x_hwp_theta,y_hwp_theta,query_hwp_wraped,sigma_hwp_theta);
    a=gauss_weighted_interp(x_hwp_cont,y_hwp_cont,query_hwp_wraped,sigma_hwp_cont);
    out_polz_state.cont.val(mask_query_no_qwp)=bound(a,0,1);
    
    
    subplot(2,3,1)
    hold on
    plot(query_hwp_wraped, out_polz_state.v.val(mask_query_no_qwp),'ko')
    hold off
    subplot(2,3,2)
    hold on
    plot(query_hwp_wraped,out_polz_state.theta.val(mask_query_no_qwp),'ko')
    hold off
    subplot(2,3,3)
    hold on
    plot(query_hwp_wraped,out_polz_state.cont.val(mask_query_no_qwp),'ko')
    hold off
    %%

    mask_query_both_with_hwp_const=pol_opts.hwp==hwp_const_val & ~isnan(pol_opts.qwp) & ~fulfilled_query;

    out_polz_state.v.val(mask_query_both_with_hwp_const)=gauss_weighted_interp(x_qwp_v,y_qwp_v,pol_opts.qwp(mask_query_both_with_hwp_const),sigma_qwp_v);
    out_polz_state.theta.val(mask_query_both_with_hwp_const)=gauss_weighted_interp(x_qwp_theta,y_qwp_theta,pol_opts.qwp(mask_query_both_with_hwp_const),sigma_qwp_theta);
    a=gauss_weighted_interp(x_qwp_cont,y_qwp_cont,pol_opts.qwp(mask_query_both_with_hwp_const),sigma_qwp_cont);
    out_polz_state.cont.val(mask_query_both_with_hwp_const)=bound(a,0,1);
    
    subplot(2,3,4)
    hold on
    plot(pol_opts.qwp(mask_query_both_with_hwp_const),  out_polz_state.v.val(mask_query_both_with_hwp_const),'ko')
    hold off
    subplot(2,3,5)
    hold on
    plot(pol_opts.qwp(mask_query_both_with_hwp_const), out_polz_state.theta.val(mask_query_both_with_hwp_const),'ko')
    hold off
    subplot(2,3,6)
    hold on
    plot(pol_opts.qwp(mask_query_both_with_hwp_const),out_polz_state.cont.val(mask_query_both_with_hwp_const),'ko')
    hold off

    %% try to see if the rest of the query points can just be extracted from the data
    if  strcmp(pol_opts.predict_method,'gauss_pref_interp')
        %mask_both_with_hwp_notconst=pol_opts.hwp~=hwp_const_val & ~isnan(pol_opts.hwp) & ~isnan(pol_opts.qwp);
        mask_unmatched=~mask_query_both_with_hwp_const & ~mask_query_no_qwp & ~fulfilled_query;
        idxs_unmatched=find(mask_unmatched);
        iimax=sum(mask_unmatched);
        for ii=1:iimax
            query_idx=idxs_unmatched(ii);
            qwp_val_query=pol_opts.qwp(query_idx);
            %get mask treating nans as equal
            is_qwp_eq= (pol_data_table.qwp_angle_deg==qwp_val_query) | (isnan(qwp_val_query) & isnan(pol_data_table.qwp_angle_deg));
            hwp_val_query=pol_opts.hwp(query_idx);
            is_hwp_eq= (pol_data_table.hwp_angle_deg==hwp_val_query) | (isnan(hwp_val_query) & isnan(pol_data_table.hwp_angle_deg));
            query_match_data= is_qwp_eq & is_hwp_eq;
            if sum(query_match_data)>0
                warning('using direct data for unmodeled points, may be missing modulo')
                data_match_idx=find(query_match_data);
                if numel(data_match_idx)>1
                    warning('multiple data matches found for query point')
                    data_match_idx=data_match_idx(1);
                end
                
                out_polz_state.cont.val(query_idx)=bound(pol_cont(data_match_idx),0,1);
                out_polz_state.theta.val(query_idx)=mod(pol_theta(data_match_idx),pi);
                out_polz_state.v.val(query_idx)=pol_v(data_match_idx);

                out_polz_state.v.unc(query_idx)=nan;
                out_polz_state.theta.unc(query_idx)=nan;
                out_polz_state.cont.unc(query_idx)=nan;
            end
        end
    end

end    
    


end
% 
% 
%     pol_data_val = [
%                        NaN    6.0000   75.0000  146.0000    0.9600   50.5000   -1.0000
%                        NaN  354.0000   89.0000  120.0000    0.9400   29.0000   -1.0000
%                        NaN  314.0000  120.0000   40.0000    0.3700  130.5000    1.0000
%                        NaN  335.0000  113.0000   84.0000    0.1300  170.0000   -1.0000
%                        NaN  344.0000  133.0000  104.0000    0.6800  190.0000   -1.0000
%                        NaN  353.0000   57.0000  298.0000    0.6500  208.0000   -1.0000
%                        NaN  347.0000   78.0000  282.0000    0.6600  194.0000   -1.0000
%                        NaN  342.0000   92.0000  278.0000    0.5200  187.0000   -1.0000
%                        NaN  340.0000   71.0000  271.0000    0.3300  181.0000   -1.0000
%                        NaN  338.5000   65.0000  266.0000    0.2300  176.5000   -1.0000
%                        NaN  338.5000   65.0000  266.0000    0.2300  176.5000   -1.0000
%                        NaN  335.0000   79.0000  264.0000    0.1600  171.0000   -1.0000
%                        NaN  332.0000   80.0000  255.0000    0.2000  165.0000    1.0000
%                        NaN  319.0000  120.0000   70.5000    0.1200  160.0000    1.0000
%                        NaN  327.0000   60.4000  248.0000    0.1100  155.0000    1.0000
%                        NaN  324.0000   87.6000  238.0000    0.2300  147.0000    1.0000
%                        NaN  320.0000  112.0000   50.0000    0.3300  140.5000    1.0000
%                        NaN  309.0000   78.0000  221.0000    0.3700  120.0000    1.0000
%                        NaN  299.0000   72.0000  188.0000    0.2500   99.5000    1.0000
%                        NaN  290.0000  107.5000  165.0000    0.0700   80.0000   -1.0000
%                        NaN  310.0000   86.0000   32.0000    0.3600  121.0000    1.0000
%                        NaN  352.0000   72.0000  112.0000    0.5400   24.5000   -1.0000
%                        NaN    3.0000   74.0000  137.0000    0.8600   46.0000   -1.0000
%                        NaN   10.0000   84.0000  154.0000    0.5200   61.0000   -1.0000
%                        NaN   26.0000   77.0000  181.0000    0.0600   92.0000    1.0000
%                        NaN    4.0000   70.1000  135.0000    0.7000  231.0000   -1.0000
%                        NaN  358.0000   78.0000  308.0000    0.9000  217.0000   -1.0000
%                        NaN  350.0000   82.0000  290.0000    0.7200  199.0000   -1.0000
%                        NaN  300.0000   72.0000  191.0000    0.1700  100.0000    1.0000
%                        NaN  320.0000   62.0000  230.0000    0.2400  140.0000    1.0000
%                        NaN  329.0000   76.0000  249.0000    0.1400  160.0000    1.0000
%                        NaN  339.0000   93.0000   89.0000    0.1900  180.0000   -1.0000
%                        NaN   15.0000  120.0000  165.0000    0.5900   70.0000   -1.0000
%                        NaN    9.5000   53.0000  150.0000    0.6400  239.5000   -1.0000
%                        NaN   14.5000   75.1000  160.0000    0.5800  249.0000   -1.0000
%                        NaN  304.0000   79.0000   20.0000    0.3400  111.0000    1.0000
%                        NaN  325.0000   98.0000   63.0000    0.1400  151.0000    1.0000
%                   280.0000  333.0000  168.0000  257.0000    0.1300  171.0000    1.0000
%                   283.0000  333.0000  195.0000   83.0000    0.6000  352.0000    1.0000
%                   283.0000  333.0000  195.0000   83.0000    0.6000  352.0000    1.0000
%                   246.0000  333.0000  127.0000   29.0000   45.0000  136.0000   -1.0000
%                   270.0000  333.0000  167.3000  248.0000    5.3700  159.0000   -1.0000
%                   286.5000  333.0000  165.0000   88.5000    2.1800  179.0000    1.0000
%                   310.0000  333.0000  140.7000  117.0000   40.8000  205.0000    1.0000
%                   286.0000  333.0000  165.0000   88.5000    2.1800  179.0000    1.0000
%                   254.0000  333.0000  115.0000  240.0000   22.1000  143.0000   -1.0000
%                   246.0000  333.0000  101.0000   50.0000   42.0000  316.0000   -1.0000
%                   234.0000  333.0000   95.0000  226.0000   82.0000  319.0000   -1.0000
%                   226.0000  333.0000   87.0000  290.0000   73.0000  199.0000   -1.0000
%                   220.0000  333.0000  110.0000  273.0000   44.1000  192.0000   -1.0000
%                   202.0000  333.0000  178.0000  262.0000   11.5000  355.0000   -1.0000
%                   187.0000  333.0000  165.0000  253.0000    0.4600  343.0000    1.0000
%                   177.0000  333.0000  150.0000  248.0000    8.7000  338.0000    1.0000
%                   162.0000  333.0000  132.0000  242.0000   35.0000  329.0000    1.0000
%                   154.0000  333.0000  113.0000  236.0000   58.0000  328.0000    1.0000
%                   130.0000  333.0000  119.0000  284.0000   31.0000   12.0000    1.0000
%                   134.0000  333.0000  156.0000  292.0000   59.0000   18.0000    1.0000
%                   138.0000  333.0000  134.0000  284.0000   58.0000    3.0000    1.0000
%                   142.0000  333.0000   77.5000  278.0000   53.0000    9.0000    1.0000
%                   146.0000  333.0000  118.0000  260.0000   74.0000  344.0000    1.0000
%                   150.0000  333.0000   90.6000  245.0000   63.9000  152.0000    1.0000
%                   260.0000  333.0000  124.0000   51.0000   18.1000  329.0000   -1.0000
%                   270.0000  340.0000  206.5000   66.0000   34.4000  156.0000   -1.0000
%                   274.0000  350.0000  205.0000  222.0000   143.200  314.0000   -1.0000
%                   268.0000  354.0000  176.0000  165.0000    87.000  252.0000   -1.0000
%                   280.0000  10.0000   148.0000  352.0000    11.200   83.0000   -1.0000
%                   290.0000   6.0000   123.5000   12.0000    49.100   281.000   -1.0000
%                   ];



% 
% 
%     
%      %% QWP,HWP,Pmax,Phi_max,p(phi_max+45),phi_max+45,Pmin,phi_min,p(phi_min+45),phi_min+45,hand
%      pol_data_table = [nan, 23,  120, 0,   62,   45,0.63,91,38.2,136 ,1   ;
%                      nan, 40,  126, 323, 60,   8,0.43,48,63,103 ,1   ;
%                      nan, 50,  115, 307, 52.4, 352,0.14,35,53,80 ,1   ;
%                      nan, 60,  100, 289, 44.8, 334,0.03,17,47,62 ,1   ;
%                      nan, 70,  120, 266, 65,   311,0.3,358,48,43 ,1   ;
%                      nan, 80,  114, 76,  44,   121,1,158,60,203 ,1   ;
%                      nan, 90,  120, 232, 55,   277,0.87,320,68,5 ,1   ;
%                      nan, 100, 121, 208, 61,   253,0.5,298,50.7,343 ,1   ;
%                      nan, 10,  100, 22,  60,   67,0.58,118,48,163 ,1   ;
%                      nan, 0,   91,  49,  61,   97,0.8,140,50.5,185 ,1   ;
%                      nan, 350, 120, 60,  72,   105,0.75,158,59,203 ,1   ;
%                      nan, 340, 112, 84,  69,   129,0.8,178,64,223 ,1   ;
%                      nan, 330, 90,  111, 43,   156,0.65,200,63,245 ,1   ;
%                      nan, 320, 117, 128, 61,   173,0.6,218,60,263 ,1   ;
%                      nan, 310, 120, 146, 60,   191,0,238,45,283 ,1   ;
%                      nan, 300, 91,  169, 53,   214,0.2,260,62,305 ,1   ;
%                      nan, 290, 120, 183, 55,   228,0.1,279,51,324 ,1   ;
%                      nan, 280, 91,  208, 46,   253,0.2,299,65,344 ,1   ;
%                      nan, 270, 113, 222, 74,   267,0.3,319,50,4 ,1   ;
%                      nan, 260, 89,  247, 49,   292,0.1,338,44.2,23 ,1   ;
%                      nan, 250, 95,  271, 38,   316,0.2,356,50,41 ,1   ;
%                      nan, 305, 100, 160, 42,   205,0.2,249,65,294 ,1   ;
%                      nan, 295, 90,  180, 46,   225,0.13,270,61,315 ,1   ;
%                      nan, 315, 116, 139, 56,   184,0.4,228,26.7,273 ,1   ;
%                      nan, 355, 88,  102, 42,   147,0.75,188,45.8,233 ,1   ;
%                      nan, 80,  110, 249, 50,   294,0.51,340,64,25 ,1   ;
%                      nan, 315, 120, 134, 68,   179,0.43,228,60,273 ,1   ;
%                      nan, 325, 120, 122, 58,   167,0.78,210,51,255 ,1   ;
%                      nan, 345, 122, 84,  40.4, 129,0.77,169,54,214 ,1   ;
%                      nan, 355, 92,  48,  61.5, 93,0.65,149,47.2,194 ,1   ;
%                      nan, 5,   91,  39,  50.1, 84,0.68,128,63,173 ,1   ;
%                      nan, 120, 90,  167, 46.2, 212,0.12,259,51,304 ,1   ;
%                      nan, 140, 85,  131, 53,   176,0.26,219,57,264 ,1   ;
%                      nan, 160, 108, 85,  60,   130,0.63,179,61,224 ,1   ;
%                      nan, 190, 112, 26,  56,   71,0.25,118,58,136 ,1   ;
%                      nan, 220, 110, 328, 54,   13,1.4,58,55,103  ,1    ;
%                      ];