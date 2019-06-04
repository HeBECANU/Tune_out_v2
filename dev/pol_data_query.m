function out_polz_state=pol_data_query(pol_opts)
% a partial rewrite of pol_data_import
%% the idea of this function is to create a bunch of differenet models of the pol state for hwp/ qwp data angles
% Inputs
% pol_opts.location = 'pre_right';%post, pre_cen, pre_left, pre_right
% pol_opts.predict = 'fit';%'interp'; %obs (obsovation) fit (pertial fit) full_fit (fit with all parameters free)
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



%%      QWP    HWP      Pmax     Phi max     Pmin     Phi min   hand
    pol_data_val = [
                       NaN    6.0000   75.0000  146.0000    0.9600   50.5000   -1.0000
                       NaN  354.0000   89.0000  120.0000    0.9400   29.0000   -1.0000
                       NaN  314.0000  120.0000   40.0000    0.3700  130.5000    1.0000
                       NaN  335.0000  113.0000   84.0000    0.1300  170.0000   -1.0000
                       NaN  344.0000  133.0000  104.0000    0.6800  190.0000   -1.0000
                       NaN  353.0000   57.0000  298.0000    0.6500  208.0000   -1.0000
                       NaN  347.0000   78.0000  282.0000    0.6600  194.0000   -1.0000
                       NaN  342.0000   92.0000  278.0000    0.5200  187.0000   -1.0000
                       NaN  340.0000   71.0000  271.0000    0.3300  181.0000   -1.0000
                       NaN  338.5000   65.0000  266.0000    0.2300  176.5000   -1.0000
                       NaN  338.5000   65.0000  266.0000    0.2300  176.5000   -1.0000
                       NaN  335.0000   79.0000  264.0000    0.1600  171.0000   -1.0000
                       NaN  332.0000   80.0000  255.0000    0.2000  165.0000    1.0000
                       NaN  319.0000  120.0000   70.5000    0.1200  160.0000    1.0000
                       NaN  327.0000   60.4000  248.0000    0.1100  155.0000    1.0000
                       NaN  324.0000   87.6000  238.0000    0.2300  147.0000    1.0000
                       NaN  320.0000  112.0000   50.0000    0.3300  140.5000    1.0000
                       NaN  309.0000   78.0000  221.0000    0.3700  120.0000    1.0000
                       NaN  299.0000   72.0000  188.0000    0.2500   99.5000    1.0000
                       NaN  290.0000  107.5000  165.0000    0.0700   80.0000   -1.0000
                       NaN  310.0000   86.0000   32.0000    0.3600  121.0000    1.0000
                       NaN  352.0000   72.0000  112.0000    0.5400   24.5000   -1.0000
                       NaN    3.0000   74.0000  137.0000    0.8600   46.0000   -1.0000
                       NaN   10.0000   84.0000  154.0000    0.5200   61.0000   -1.0000
                       NaN   26.0000   77.0000  181.0000    0.0600   92.0000    1.0000
                       NaN    4.0000   70.1000  135.0000    0.7000  231.0000   -1.0000
                       NaN  358.0000   78.0000  308.0000    0.9000  217.0000   -1.0000
                       NaN  350.0000   82.0000  290.0000    0.7200  199.0000   -1.0000
                       NaN  300.0000   72.0000  191.0000    0.1700  100.0000    1.0000
                       NaN  320.0000   62.0000  230.0000    0.2400  140.0000    1.0000
                       NaN  329.0000   76.0000  249.0000    0.1400  160.0000    1.0000
                       NaN  339.0000   93.0000   89.0000    0.1900  180.0000   -1.0000
                       NaN   15.0000  120.0000  165.0000    0.5900   70.0000   -1.0000
                       NaN    9.5000   53.0000  150.0000    0.6400  239.5000   -1.0000
                       NaN   14.5000   75.1000  160.0000    0.5800  249.0000   -1.0000
                       NaN  304.0000   79.0000   20.0000    0.3400  111.0000    1.0000
                       NaN  325.0000   98.0000   63.0000    0.1400  151.0000    1.0000
                  280.0000  333.0000  168.0000  257.0000    0.1300  171.0000    1.0000
                  283.0000  333.0000  195.0000   83.0000    0.6000  352.0000    1.0000
                  283.0000  333.0000  195.0000   83.0000    0.6000  352.0000    1.0000
                  246.0000  333.0000  127.0000   29.0000   45.0000  136.0000   -1.0000
                  270.0000  333.0000  167.3000  248.0000    5.3700  159.0000   -1.0000
                  286.5000  333.0000  165.0000   88.5000    2.1800  179.0000    1.0000
                  310.0000  333.0000  140.7000  117.0000   40.8000  205.0000    1.0000
                  286.0000  333.0000  165.0000   88.5000    2.1800  179.0000    1.0000
                  254.0000  333.0000  115.0000  240.0000   22.1000  143.0000   -1.0000
                  246.0000  333.0000  101.0000   50.0000   42.0000  316.0000   -1.0000
                  234.0000  333.0000   95.0000  226.0000   82.0000  319.0000   -1.0000
                  226.0000  333.0000   87.0000  290.0000   73.0000  199.0000   -1.0000
                  220.0000  333.0000  110.0000  273.0000   44.1000  192.0000   -1.0000
                  202.0000  333.0000  178.0000  262.0000   11.5000  355.0000   -1.0000
                  187.0000  333.0000  165.0000  253.0000    0.4600  343.0000    1.0000
                  177.0000  333.0000  150.0000  248.0000    8.7000  338.0000    1.0000
                  162.0000  333.0000  132.0000  242.0000   35.0000  329.0000    1.0000
                  154.0000  333.0000  113.0000  236.0000   58.0000  328.0000    1.0000
                  130.0000  333.0000  119.0000  284.0000   31.0000   12.0000    1.0000
                  134.0000  333.0000  156.0000  292.0000   59.0000   18.0000    1.0000
                  138.0000  333.0000  134.0000  284.0000   58.0000    3.0000    1.0000
                  142.0000  333.0000   77.5000  278.0000   53.0000    9.0000    1.0000
                  146.0000  333.0000  118.0000  260.0000   74.0000  344.0000    1.0000
                  150.0000  333.0000   90.6000  245.0000   63.9000  152.0000    1.0000
                  260.0000  333.0000  124.0000   51.0000   18.1000  329.0000   -1.0000
                  270.0000  340.0000  206.5000   66.0000   34.4000  156.0000   -1.0000
                  274.0000  350.0000  205.0000  222.0000   143.200  314.0000   -1.0000
                  268.0000  354.0000  176.0000  165.0000    87.000  252.0000   -1.0000
                  280.0000  10.0000   148.0000  352.0000    11.200   83.0000   -1.0000
                  290.0000   6.0000   123.5000   12.0000    49.100   281.000   -1.0000
                  ];
    %if we want to use the observation method
    pol_v = pol_data_val(:,7).*2.*sqrt(pol_data_val(:,3).*pol_data_val(:,5))...
        ./(pol_data_val(:,3)+pol_data_val(:,5));%the V parameter for each run
    pol_cont = (pol_data_val(:,3)-pol_data_val(:,5))./(pol_data_val(:,3)+pol_data_val(:,5));%contrast
    pol_theta = pol_data_val(:,6).*pi/180; %using min pow angle

    
if strcmp(pol_opts.predict,'full_fit') 
    % make a model for the case of having qwp+hwp, and just qwp
    
    idx_qwp_and_hwp=~isnan(pol_data_val(:,1)) & ~isnan(pol_data_val(:,2));
    hwp_const_val= mode(pol_data_val(~isnan(pol_data_val(:,1)),2))
    idx_qwp_and_hwp_const=~isnan(pol_data_val(:,1)) & pol_data_val(:,2)==hwp_const_val;
    idx_hwp=~isnan(pol_data_val(:,2)) & isnan(pol_data_val(:,1));
    idx_qwp=isnan(pol_data_val(:,2)) & ~isnan(pol_data_val(:,1));
    if sum(idx_qwp)>0
        error('not set up to deal with just qwp')
    end
    stfig('polz fits');
    clf
    
    %% fit the 4th stokes parameter for the hwp data
    beta0 = [-0.18297,1.5186,-0.033813];
    fit_V_hwp_only = fitnlm(pol_data_val(idx_hwp,2),pol_v(idx_hwp),sin_mdl_fixed_freq4,beta0);
    subplot(2,3,1)
    plot(pol_data_val(idx_hwp,2),pol_v(idx_hwp),'x')
    hold on
    hwp_samp_vec=col_vec(linspace(min(pol_data_val(idx_hwp,2)),max(pol_data_val(idx_hwp,2)),1e3));
    [v_samp_fit_val,V_samp_fit_ci]=predict(fit_V_hwp_only,hwp_samp_vec);
    plot(hwp_samp_vec,v_samp_fit_val,'-k')
    plot(hwp_samp_vec,V_samp_fit_ci,'-b')
    plot(hwp_samp_vec,sin_mdl_fixed_freq4(beta0,hwp_samp_vec),'-r')
    hold off
    
        
    %fit_theta_hwp = fitnlm(mod(hwp_ang+21,90),mod(theta(1:37),pi),lin_mdl,beta0);
    %beta0 = [1.0,2.0139,-0.41317,-0.023839];
    %fit_d_p_qwp = fitnlm(qwp_ang,d_p(38:63),sin_mdl_abs,beta0);
    
    %% fit the theta 2nd,3rd stokes parameter angle for the hwp data
    hwp_wraped=mod(pol_data_val(idx_hwp,2),90);
    theta_wraped=mod(pol_theta(idx_hwp),pi);
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
    
    %% fit the theta 2nd,3rd stokes parameter contrast for the hwp data

    poly_mdl=@(b,x) taylor_series(x,b,90/2);
    %beta0 = [1,0,-0.001,0.001,1e-6,1e-6,1e-7,1e-8];
    beta0=[9.923421e-01,3.107072e-04,5.344973e-05,-8.093548e-06,-1.188569e-06,1.759216e-07,1.218760e-08,-2.301748e-09];
     opt = statset('TolFun',1e-10,'TolX',1e-10,...
            'MaxIter',1e4,... %1e4
            'UseParallel',1);
    %coef_names={'amp','phase','offset'};
    fit_cont_hwp_only = fitnlm(hwp_wraped,pol_cont(idx_hwp),poly_mdl,beta0,...
        'options',opt);         %'CoefficientNames',coef_names,...
    subplot(2,3,3)
    plot(hwp_wraped,pol_cont(idx_hwp),'x')
    hold on
    hwp_samp_vec=col_vec(linspace(min(hwp_wraped),max(hwp_wraped),1e3));
    [theta_samp_fit_val,theta_samp_fit_ci]=predict(fit_cont_hwp_only,hwp_samp_vec);
    plot(hwp_samp_vec,theta_samp_fit_val,'-k')
    plot(hwp_samp_vec,theta_samp_fit_ci,'-b')
    %plot(hwp_samp_vec,poly_mdl(beta0,hwp_samp_vec),'-r')
    hold off
    
    
    %% fit the 4th stokes parameter for the qwp data (with hwp=333)
     
    % fit only the dominant hwp angle here
    %idx_qwp_and_hwp_333
    beta0 = [-1.0002,-0.41317,-0.023839];
    fit_V_qwp_and_hwp_const = fitnlm(pol_data_val(idx_qwp_and_hwp_const,1),pol_v(idx_qwp_and_hwp_const),sin_mdl_fixed_freq2,beta0);
    subplot(2,3,4)
    plot(pol_data_val(idx_qwp_and_hwp_const,1),pol_v(idx_qwp_and_hwp_const),'x')
    hold on
    qwp_angle_samp_fit_hwp_const=col_vec(linspace(min(pol_data_val(idx_qwp_and_hwp_const,1)),max(pol_data_val(idx_qwp_and_hwp_const,1)),1e3));
    [v_samp_fit_val,V_samp_fit_ci]=predict(fit_V_qwp_and_hwp_const,qwp_angle_samp_fit_hwp_const);
    plot(qwp_angle_samp_fit_hwp_const,v_samp_fit_val,'-k')
    plot(qwp_angle_samp_fit_hwp_const,V_samp_fit_ci,'-b')
    plot(qwp_angle_samp_fit_hwp_const,sin_mdl_fixed_freq2(beta0,qwp_angle_samp_fit_hwp_const),'-r')
    hold off
    xlabel('qwp angle')
    ylabel('4th stokes')
    
    %% fit the theta 2nd,3rd stokes angle parameter for the qwp data (with hwp=333)
    qwp_wraped=mod(pol_data_val(idx_qwp_and_hwp_const,1),90);
    theta_wraped=mod(pol_theta(idx_qwp_and_hwp_const),pi);
    [~,sort_idx]=sort(qwp_wraped);
    theta_wraped(sort_idx)=unwrap(theta_wraped(sort_idx)*8)/8;
    
    beta0=[0.5,-1,2.5];
    opt = statset('TolFun',1e-10,'TolX',1e-10,...
            'MaxIter',1e4,... %1e4
            'UseParallel',1);
    %sin_mdl = @(b,x) b(1).*sin(b(2).*x.*pi/180+b(3))+b(4);
    coef_names={'amp','phase','offset'};
    fit_theta_qwp_and_hwp_const = fitnlm(pol_data_val(idx_qwp_and_hwp_const,1),theta_wraped,sin_mdl_fixed_freq4,beta0,...
        'CoefficientNames',coef_names,...
        'options',opt);
    subplot(2,3,5)
    plot(pol_data_val(idx_qwp_and_hwp_const,1),theta_wraped,'x')
    hold on
    v_hwp_angle_samp_fit_hwp_const=col_vec(linspace(min(pol_data_val(idx_qwp_and_hwp_const,1)),max(pol_data_val(idx_qwp_and_hwp_const,1)),1e3));
    [v_samp_fit_val,V_samp_fit_ci]=predict(fit_theta_qwp_and_hwp_const,v_hwp_angle_samp_fit_hwp_const);
    plot(v_hwp_angle_samp_fit_hwp_const,v_samp_fit_val,'-k')
    plot(v_hwp_angle_samp_fit_hwp_const,V_samp_fit_ci,'-b')
    plot(v_hwp_angle_samp_fit_hwp_const,sin_mdl_fixed_freq4(beta0,v_hwp_angle_samp_fit_hwp_const),'-r')
    hold off
    title('theta')
    
    %% fit the theta 2nd,3rd stokes parameter contrast for the qwp data (with hwp=333)
    
    beta0 = [0.3,1.5186,0.5];
     opt = statset('TolFun',1e-10,'TolX',1e-10,...
            'MaxIter',1e4,... %1e4
            'UseParallel',1);
    coef_names={'amp','phase','offset'};
    fit_cont_qwp_const_hwp = fitnlm(pol_data_val(idx_qwp_and_hwp_const,1),pol_cont(idx_qwp_and_hwp_const),sin_mdl_fixed_freq4,beta0,...
        'CoefficientNames',coef_names,...
        'options',opt);
    subplot(2,3,6)
    plot(pol_data_val(idx_qwp_and_hwp_const,1),pol_cont(idx_qwp_and_hwp_const),'x')
    hold on
    hwp_samp_vec=col_vec(linspace(min(pol_data_val(idx_qwp_and_hwp_const,1)),max(pol_data_val(idx_qwp_and_hwp_const,1)),1e3));
    [theta_samp_fit_val,theta_samp_fit_ci]=predict(fit_cont_qwp_const_hwp,hwp_samp_vec);
    plot(hwp_samp_vec,theta_samp_fit_val,'-k')
    plot(hwp_samp_vec,theta_samp_fit_ci,'-b')
    plot(hwp_samp_vec,sin_mdl_fixed_freq4(beta0,hwp_samp_vec),'-r')
    hold off
    
    %% now calulate at the query points
    only_qwp=isnan(pol_opts.hwp) & ~isnan( pol_opts.qwp);
    if sum( only_qwp)>0
        warning('not set up to query quater wp by itself')
        out_theta_v(only_qwp,1)=nan;
        out_theta_v(only_qwp,2)=nan;
    end
    
    
    idx_no_qwp=~isnan(pol_opts.hwp) & isnan(pol_opts.qwp);
    query_hwp_wraped=mod(pol_opts.hwp(idx_no_qwp),90);
    
    
    [a,b]=predict(fit_theta_hwp_only,pol_opts.hwp(idx_no_qwp),'Alpha',1-erf(1/sqrt(2)));
    out_polz_state.theta.val(idx_no_qwp)=a;
    out_polz_state.theta.unc(idx_no_qwp)=range(b,2)/2;
    
    
    [a,b]=predict(fit_V_hwp_only,query_hwp_wraped,'Alpha',1-erf(1/sqrt(2)));
    out_polz_state.v.val(idx_no_qwp)=a;
    out_polz_state.v.unc(idx_no_qwp)=range(b,2)/2;
    
    [a,b]=predict(fit_cont_hwp_only,query_hwp_wraped,'Alpha',1-erf(1/sqrt(2)));
    out_polz_state.cont.val(idx_no_qwp)=a;
    out_polz_state.cont.unc(idx_no_qwp)=range(b,2)/2;
    
    subplot(2,3,1)
    hold on
    plot(pol_opts.hwp(idx_no_qwp),out_polz_state.theta.val(idx_no_qwp),'ko')
    hold off
    subplot(2,3,2)
    hold on
    plot(query_hwp_wraped, out_polz_state.v.val(idx_no_qwp),'ko')
    hold off
    subplot(2,3,3)
    hold on
    plot(query_hwp_wraped,out_polz_state.cont.val(idx_no_qwp),'ko')
    hold off
    %%


    idx_both_with_hwp_const=pol_opts.hwp==hwp_const_val & ~isnan(pol_opts.qwp);


    [a,b]=predict(fit_V_qwp_and_hwp_const,pol_opts.qwp(idx_both_with_hwp_const),'Alpha',1-erf(1/sqrt(2)));
    out_polz_state.v.val(idx_both_with_hwp_const)=a;
    out_polz_state.v.unc(idx_both_with_hwp_const)=range(b,2)/2;
    [a,b]=predict(fit_theta_qwp_and_hwp_const,pol_opts.qwp(idx_both_with_hwp_const),'Alpha',1-erf(1/sqrt(2)));
    out_polz_state.theta.val(idx_both_with_hwp_const)=a;
    out_polz_state.theta.unc(idx_both_with_hwp_const)=range(b,2)/2;
    
    [a,b]=predict(fit_cont_qwp_const_hwp,pol_opts.qwp(idx_both_with_hwp_const),'Alpha',1-erf(1/sqrt(2)));
    out_polz_state.cont.val(idx_both_with_hwp_const)=a;
    out_polz_state.cont.unc(idx_both_with_hwp_const)=range(b,2)/2;
    
    subplot(2,3,4)
    hold on
    plot(pol_opts.qwp(idx_both_with_hwp_const),  out_polz_state.v.val(idx_both_with_hwp_const),'ko')
    hold off
    subplot(2,3,5)
    hold on
    plot(pol_opts.qwp(idx_both_with_hwp_const), out_polz_state.theta.val(idx_both_with_hwp_const),'ko')
    hold off
    subplot(2,3,6)
    hold on
    plot(pol_opts.qwp(idx_both_with_hwp_const),out_polz_state.cont.val(idx_both_with_hwp_const),'ko')
    hold off

    
    idx_both_with_hwp_notconst=pol_opts.hwp~=hwp_const_val & ~isnan(pol_opts.hwp) & ~isnan(pol_opts.qwp);
    iimax=sum(idx_both_with_hwp_notconst);
    if iimax>1
        error('havent done this yet')
    end
    
end






end