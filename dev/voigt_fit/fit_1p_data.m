function out_st=fit_1p_data(opts)
predictor=opts.predictor;
response=opts.response;

if ~isfield(opts,'do_plots')
    opts.do_plots=true;
end

if ~isfield(opts,'sigma_outlier')
    opts.sigma_outlier=4;
end

provided_ste=false;
if isfield(opts,'response_err')
    provided_ste=true;
end

fit_noise=false;
if isfield(opts,'response_noise')
    fit_noise=true;
end


out_st=[];
%% fit to a lorentzian to start with
%ini_fit_fun=@(b,x) b(3)*exp(-(1/2)*((x-b(2))./b(1)).^2)+b(4);
ini_fit_fun=@(b,x) lorentzian_function_1d(x,b(1),b(2),b(3),b(4))+b(5)*x+b(6)*x.^2;
cof_names={'gamma','mu','amp','offset','grad','curve'};


% do a robust fit to get the global
gf_opt=[];
xscale=range(predictor);
xmean=mean(predictor);
yscale=range(response);
ymean=mean(response);
gf_opt.domain=[[0.02,0.5]*xscale;...   %sigma
               [-0.3,0.3]*xscale+xmean;...       %mu
               [-1,1]*yscale;...     %amp
               [-1,1]*yscale+ymean;... %offset
               [-1,1]*2*yscale/xscale;... %grad
               [-1,1]*2*(yscale^2)*xscale;... %offset
               ];        
gf_opt.start=[(1/2)*range(predictor),mean(predictor),range(response),mean(response),0,0];
gf_opt.rmse_thresh=0.1*yscale;
gf_opt.plot=false;
gf_opt.level=2;
gf_out=global_fit(predictor,response,ini_fit_fun,gf_opt);

% now do a normal fit using the global fit as a start pt
beta0=gf_out.params;
nlmopt = statset('UseParallel',1);
% 'TolFun',1e-10,'TolX',1e-10,...
%     'MaxIter',1e4,... %1e4
%     );
% weights=ones(size(predictor));
% %'Weights',weights
fit_obj_lorentzian=fitnlm(predictor,response,ini_fit_fun,beta0,...
    'options',nlmopt,...
    'CoefficientNames',cof_names);

% %%
% % lets split up the fwhm that was just found into the gauss and lorentzian component
% % ill bastadrdize the forumula from https://doi.org/10.1016%2F0022-4073%2877%2990161-3
% % implemented in voigt_approx_fwhm
% % 0.5346*fwhm_l+sqrt(0.2166*(fwhm_l^2)+fwhm_g^2)
% gamma_ini_fit=fit_obj_lorentzian.Coefficients.Estimate(1);
% fwhm_g=2*gamma_ini_fit*sqrt(2*log(2));
% % now i split up the fwhm into the gauss, and lorentzian componets of the voit
% % i have arbitrarily chosen to split them so they each contribute half to the fwhm
% % this is a pretty good place to start
% fwhm_v_g=fwhm_g*0.8;
% fwhm_v_l= 0.4*fwhm_g; 
% gamma_start=fwhm_v_l/2;
% sigma_start=fwhm_v_g/(sqrt(2*log(2))*2);
% mu_start=fit_obj_lorentzian.Coefficients.Estimate(2);
% amp_start=fit_obj_lorentzian.Coefficients.Estimate(3);
% offset_start=fit_obj_lorentzian.Coefficients.Estimate(4);
% % diagnose how the splitting went
% % fwhm_g
% % voigt_approx_fwhm(sigma_start,gamma_start)
% %%
% 
% 
% vo_fit_fun = @(b,x) voigt_function_1d(x,b(1),b(2),b(3),b(4),b(5),'norm','amp','method','approx')+b(6)*x+b(7)*x.^2;
% cof_names={'sigma','gamma','mu','amp','offset','grad','curve'};
% %vo_fit_fun = @(b,x) voigt_function_1d(x,b(1),b(2),b(3),b(4),b(5),'norm','amp','method','fadd');
% 
%   
% beta0=[sigma_start,gamma_start,mu_start,amp_start,offset_start,0,0];
% 
% nlmopt = statset('UseParallel',1);
% % 'TolFun',1e-10,'TolX',1e-10,...
% %     'MaxIter',1e4,... %1e4
% %     );
% % weights=ones(size(predictor));
% % %'Weights',weights
% 
% fitobj=fitnlm(predictor,response,vo_fit_fun,beta0,...
%     'options',nlmopt,...
%     'CoefficientNames',cof_names,'ErrorModel','combined');
% 
%%
fitobj=fit_obj_lorentzian;



%% look for outliers
out_st.fit_no_cull.voigt=fitobj;
out_st.fit_no_cull.lorentzian=fit_obj_lorentzian;

if ~isinf(opts.sigma_outlier)
    [~,yci_cull_lim]=fitobj.predict(predictor,'Alpha',1-erf(opts.sigma_outlier/sqrt(2)),'Prediction','observation'); %
    is_in_ci_mask=response>yci_cull_lim(:,1) & response<yci_cull_lim(:,2);
else
    is_in_ci_mask=true(size(response));
end



if opts.do_plots
    if any(~is_in_ci_mask)
        stfig('fit results with outliers');
    else
        stfig('fit results');
    end
    clf
    subplot(10,1,[1,7])
    hold on
    %fitobj=fit_obj_ini;

    y_offset=fitobj.Coefficients{'offset','Estimate'};
    x_offset=fitobj.Coefficients{'mu','Estimate'};
    
    xplotvalues=linspace(min(predictor),max(predictor),1e4);
    xplotvalues=col_vec(xplotvalues);
    %[prediction,ci]=predict(fitobject,tplotvalues); %'Alpha',1-erf(1/sqrt(2)),'Prediction','observation'
    [amp_pred,ci_obs]=fitobj.predict(xplotvalues,'Alpha',1-erf(1/sqrt(2)),'Prediction','observation'); %
    [~,ci_curve]=fitobj.predict(xplotvalues,'Alpha',1-erf(1/sqrt(2)),'Prediction','curve'); %

    %prediction=prediction(1:end-1);
    %

    
    plot_ci_obs=true;
    shaded_lines=true;
    if plot_ci_obs
        if shaded_lines
            color_shaded=[0.9,0.9,0.9];
            patch([(xplotvalues-x_offset)', fliplr((xplotvalues-x_offset)')],...
                [ci_obs(:,1)', fliplr(ci_obs(:,2)')]-y_offset, color_shaded,'EdgeColor','none')  %[1,1,1]*0.80
            hold on
        else
            color_shaded=[0,1,0];
            plot(xplotvalues-x_offset,ci_obs(:,1)-y_offset,'-','LineWidth',1.5,'Color',color_shaded)
            hold on
            plot(xplotvalues-x_offset,ci_obs(:,2)-y_offset,'-','LineWidth',1.5,'Color',color_shaded)
        end  
    end
    
    plot_ci_curve=true;
    shaded_lines=true;
    if plot_ci_curve
        if shaded_lines
            color_shaded=[0.9,1,1];
            patch([(xplotvalues-x_offset)', fliplr((xplotvalues-x_offset)')],...
                [ci_curve(:,1)', fliplr(ci_curve(:,2)')]-y_offset, color_shaded,'EdgeColor','none')  %[1,1,1]*0.80
            hold on
        else
            color_shaded=[0,1,0];
            plot(xplotvalues-x_offset,ci_curve(:,1)-y_offset,'-','LineWidth',1.5,'Color',color_shaded)
            hold on
            plot(xplotvalues-x_offset,ci_curve(:,2)-y_offset,'-','LineWidth',1.5,'Color',color_shaded)
        end  
    end

    plot(xplotvalues-x_offset,amp_pred-y_offset,'k-','LineWidth',1.0)
    if provided_ste
        errorbar(predictor-x_offset,response-y_offset,opts.response_err,opts.response_err,'ko'...
        ,'CapSize',0,'MarkerSize',3,...
        'LineWidth',1.5,...
        'MarkerEdgeColor',[0.3,0.3,0.8],...
        'MarkerFaceColor',[0.45,0.45,1])
    end
     %plot(predictor,response-y_offset,'k.',...
     %    'MarkerEdgeColor','red',... 
     %   )
    
    plot(predictor(~is_in_ci_mask)-x_offset,response(~is_in_ci_mask)-y_offset,'xr')
    hold off 
    
    ylabel('PD Signal (V)')
    xmin_max=min_max_vec(predictor-x_offset);
    xrange=range(xmin_max);
    xticks([(ceil(xmin_max(1)/10)*10):10:0,10:10:(floor(xmin_max(2)/10)*10)])
    set(gca,'Xticklabel',[]) 
    %xlabel('Frequency (MHz)')
    xlim(xmin_max+[-0.2,0.2])
    xl=xlim;
    xtic=xticks;
    yl=ylim;
    ylim([-0.2,yl(2)])
    set(gca,'fontsize', 12)
    box on
    set(gca,'LineWidth',1)
    set(gca,'XColor','k')
    set(gca,'YColor','k')
    
    %% plot residuals
    subplot(10,1,[8,10])
    plot_ci_obs=true;
    shaded_lines=true;
    if plot_ci_obs
        if shaded_lines
            color_shaded=[0.9,0.9,0.9];
            patch([(xplotvalues-x_offset)', fliplr((xplotvalues-x_offset)')],...
                [(ci_obs(:,1)-amp_pred)', fliplr((ci_obs(:,2)-amp_pred)')],...
                color_shaded,'EdgeColor','none')  %[1,1,1]*0.80
            hold on
        else
            color_shaded=[0,1,0];
            plot(xplotvalues-x_offset,ci_obs(:,1)-amp_pred,'-','LineWidth',1.5,'Color',color_shaded)
            hold on
            plot(xplotvalues-x_offset,ci_obs(:,2)-amp_pred,'-','LineWidth',1.5,'Color',color_shaded)
        end  
    end
    
    plot_ci_curve=true;
    shaded_lines=false;
    if plot_ci_curve
        if shaded_lines
            color_shaded=[0.9,1,1];
            patch([(xplotvalues-x_offset)', fliplr((xplotvalues-x_offset)')],...
                [(ci_curve(:,1)-amp_pred)', fliplr((ci_curve(:,2)-amp_pred)')],...
                color_shaded,'EdgeColor','none')  %[1,1,1]*0.80
            hold on
        else
            color_shaded=[0,0.4,0];
            plot(xplotvalues-x_offset,ci_curve(:,1)-amp_pred,':','LineWidth',1.5,'Color',color_shaded)
            hold on
            plot(xplotvalues-x_offset,ci_curve(:,2)-amp_pred,':','LineWidth',1.5,'Color',color_shaded)
        end  
    end
    
    plot(xplotvalues-x_offset,amp_pred*0,'k-','LineWidth',1.0)
    % calculate the predicted values at predictor points
    amp_pred=fitobj.predict(predictor); %
    if provided_ste
        ebar=errorbar(predictor-x_offset,response-amp_pred,opts.response_err,opts.response_err,'ko'...
        ,'CapSize',0,'MarkerSize',3,...
        'LineWidth',1.5,...
        'MarkerEdgeColor',[0.3,0.3,0.8],...
        'MarkerFaceColor',[0.45,0.45,1],...
        'Color',[0.45,0.45,1]);
    end
     %plot(predictor,response-y_offset,'k.',...
     %    'MarkerEdgeColor','red',... 
     %   )
    
    plot(predictor(~is_in_ci_mask)-x_offset,response(~is_in_ci_mask)-amp_pred(~is_in_ci_mask),'xr')
    hold off 
    xticks(xtic)
    xlim(xl)
    ylabel('Fit Residual (V)')
    xlabel('Frequency (MHz)')
    set(gca,'fontsize', 12)
    pause(1e-9)
    box on
    set(gca,'LineWidth',1)
    set(gca,'XColor','k')
    set(gca,'YColor','k')
   
    
end

%% fit noise
if fit_noise
    %      % I use a fit to the noise to procude another measurment of the offset & the noise in the
    %         % laser
    %         % use all the previously fitted params(such as the width) except for the offset
    %         %then use the derivative of this as the laser frequency noise componet of the PMT noise
    %         % The noise model we use is laser_freq_noise*d signal/d freq + gain(volts/photon) * sqrt(signal)+
    %         % offset
    %         % this is really usefull because it gives an estimate of the system gain and the laser noise
    %         % in CCD this is called the photon transfer technique
    %         %http://hosting.astro.cornell.edu/academics/courses/astro3310/Books/Janesick_PhotonTransfer_SPIE1987.pdf
    %         %http://spiff.rit.edu/classes/phys445/lectures/gain/gain.html
    %         % derviing the noise sources
    %         % sigma Signal = sqrt(sigma Background^2 + sigma Probe^2 )
    %         % Signal =Background + Probe
    %         % Probe= Signal - Background
    %         % sigma Signal^2 = sigma Background^2 + sigma Probe^2
    %         % from shot noise sigma Probe =sqrt (Probe)/ sqrt(k) (gain in photons/(volt·s)
    %         % sigma Signal^2 = sigma Background^2 + (sqrt (Probe)/ k )^2
    %         % sigma Signal = sqrt((sqrt (Signal-Background)/ k )^2 + sigma Background^2)
    %         % sigma Signal = sqrt(( (Signal-Background)/ k ) + sigma Background^2)




    % sqrt(amp_fit_fun(x,b(2)))./sqrt(abs(b(3))) +b(4)   
%     noise_fit_fun = @(b,x) b(1) *abs(lorentzian_function_1d(x,fitobj.Coefficients{'gamma','Estimate'},...
%                                     fitobj.Coefficients{'mu','Estimate'},'norm','amp','derivative',1)) ...
%                             +abs(b(2))*sqrt(voigt_function_1d(x,fitobj.Coefficients{'sigma','Estimate'},...
%                                     fitobj.Coefficients{'gamma','Estimate'},fitobj.Coefficients{'mu','Estimate'},...
%                                     fitobj.Coefficients{'amp','Estimate'},fitobj.Coefficients{'sigma','Estimate'},...
%                                     'norm','amp','method','approx'))...
%                              + b(3);
    % we have lots of components
    % laser intensity noise  propto signal
    % signal shot noide      propto sqrt signal
    % laser freq noise      propto abs deriv signal
    % background noise
                     

%     noise_fit_fun = @(b,x) sqrt( (b(1) *abs(voigt_function_1d(x,fitobj.Coefficients{'sigma','Estimate'},fitobj.Coefficients{'gamma','Estimate'},...
%                                 fitobj.Coefficients{'mu','Estimate'},fitobj.Coefficients{'amp','Estimate'},'norm','amp','derivative',1))).^2 ...
%                         +(abs(b(2))*sqrt(voigt_function_1d(x,fitobj.Coefficients{'sigma','Estimate'},...
%                                 fitobj.Coefficients{'gamma','Estimate'},fitobj.Coefficients{'mu','Estimate'},...
%                                 fitobj.Coefficients{'amp','Estimate'},fitobj.Coefficients{'sigma','Estimate'},...
%                                 'norm','amp','method','approx'))).^2 ...
%                         + (abs(b(3))*voigt_function_1d(x,fitobj.Coefficients{'sigma','Estimate'},...
%                                 fitobj.Coefficients{'gamma','Estimate'},fitobj.Coefficients{'mu','Estimate'},...
%                                 fitobj.Coefficients{'amp','Estimate'},fitobj.Coefficients{'sigma','Estimate'},...
%                                 'norm','amp','method','approx')).^2 ...
%                          + b(4)^2 );                        
noise_fit_wrap=@(b,x) noise_fit_fun (b,x,fitobj);
    cof_names={'freq_n','amp_shot_n','back_n'};                                      
                         
    response_noise=opts.response_noise;
    beta0=[0.1,0.2,mean(response_noise)];

    nlmopt = statset('UseParallel',1,...
        'MaxIter',1e3);
    % 'TolFun',1e-10,'TolX',1e-10,...
    %     'MaxIter',1e4,... %1e4
    %     );
    % weights=ones(size(predictor));
    % %'Weights',weights

    fitobj_noise=fitnlm(predictor,response_noise,noise_fit_wrap,beta0,...
        'options',nlmopt,...
        'CoefficientNames',cof_names); %'ErrorModel','combined'
    out_st.fit_noise_no_cull=fitobj_noise;  
    if opts.do_plots
        stfig('noise fit results');
        clf
        xplotvalues=linspace(min(predictor),max(predictor),1e4);
        xplotvalues=col_vec(xplotvalues);
        y_offset=0
        %[prediction,ci]=predict(fitobject,tplotvalues); %'Alpha',1-erf(1/sqrt(2)),'Prediction','observation'
        [amp_pred,ci_obs]=fitobj_noise.predict(xplotvalues,'Alpha',1-erf(1/sqrt(2)),'Prediction','observation'); %
        [~,ci_curve]=fitobj_noise.predict(xplotvalues,'Alpha',1-erf(1/sqrt(2)),'Prediction','curve'); %
        plot_ci_obs=true;
        shaded_lines=true;
        if plot_ci_obs
            if shaded_lines
                color_shaded=[0.9,1,1];
                patch([(xplotvalues)', fliplr((xplotvalues-x_offset)')],...
                    [ci_obs(:,1)', fliplr(ci_obs(:,2)')]-y_offset, color_shaded,'EdgeColor','none')  %[1,1,1]*0.80
                hold on
            else
                color_shaded=[0,1,0];
                plot(xplotvalues-x_offset,ci_obs(:,1)-y_offset,'-','LineWidth',1.5,'Color',color_shaded)
                hold on
                plot(xplotvalues-x_offset,ci_obs(:,2)-y_offset,'-','LineWidth',1.5,'Color',color_shaded)
            end  
        end
    
        plot_ci_curve=false;
        shaded_lines=false;
        if plot_ci_curve
            if shaded_lines
                color_shaded=[0.9,1,1];
                patch([(xplotvalues-x_offset)', fliplr((xplotvalues-x_offset)')],...
                    [ci_curve(:,1)', fliplr(ci_curve(:,2)')]-y_offset, color_shaded,'EdgeColor','none')  %[1,1,1]*0.80
                hold on
            else
                color_shaded=[0,1,0];
                plot(xplotvalues-x_offset,ci_curve(:,1)-y_offset,'-','LineWidth',1.5,'Color',color_shaded)
                hold on
                plot(xplotvalues-x_offset,ci_curve(:,2)-y_offset,'-','LineWidth',1.5,'Color',color_shaded)
            end  
        end

        plot(xplotvalues-x_offset,amp_pred-y_offset,'k-','LineWidth',1.0)
        % http://davidmlane.com/hyperstat/A19196.html
        % standard err of the standard deviation
        response_noise_err=0.17*response_noise/sqrt(opts.num_samples_per_pt);
        if provided_ste
            errorbar(predictor-x_offset,response_noise-y_offset,response_noise_err,response_noise_err,'ko'...
            ,'CapSize',0,'MarkerSize',3,...
            'LineWidth',1.5,...
            'MarkerEdgeColor',[0.3,0.3,0.8],...
            'MarkerFaceColor',[0.45,0.45,1])
        end
         %plot(predictor,response-y_offset,'k.',...
         %    'MarkerEdgeColor','red',... 
         %   )

        plot(predictor(~is_in_ci_mask)-x_offset,response(~is_in_ci_mask)-y_offset,'xr')
        hold off 
    
        ylabel('PMT Current noise (nA)')
        xmin_max=min_max_vec(predictor-x_offset);
        xrange=range(xmin_max);
        xticks([(ceil(xmin_max(1)/2)*2):2:0,2:2:(2*floor(xmin_max(2)/2))])
        %set(gca,'Xticklabel',[]) 
        %xlabel('Frequency (MHz)')
        xlim(xmin_max+[-0.2,0.2])
        xl=xlim;
        xtic=xticks;
        yl=ylim;
        %ylim([-0.2,yl(2)])
        set(gca,'fontsize', 12)
        box on
        set(gca,'LineWidth',1)
        set(gca,'XColor','k')
        set(gca,'YColor','k')
        xlabel('Frequency (MHz)')
    end
end
%%


if any(~is_in_ci_mask)
    out_st.fit_no_cull.obj=fitobj;
    opts.predictor=opts.predictor(is_in_ci_mask);
    opts.response=opts.response(is_in_ci_mask);
    if fit_noise
        opts.response_noise=opts.response_noise(is_in_ci_mask);
    end
    opts.response_err=opts.response_err(is_in_ci_mask);
    if numel(opts.predictor)>20
       % run the fit again on the culled data (without culling)
       opts.sigma_outlier=inf;
       out_st_culled=fit_2p_data_with_a_voigt(opts);
       out_st.fit_cull=out_st_culled.fit_no_cull;
       if fit_noise
        out_st.fit_noise_cull=out_st_culled.fit_noise_no_cull;
       end
    else
       out_st.fit_cull=[];
    end
else
    out_st.fit_cull=out_st.fit_no_cull;
    if fit_noise
        out_st.fit_noise_cull=out_st.fit_noise_no_cull;
    end
end
out_st.data_cull=[];
out_st.data_cull.predictor=opts.predictor;
out_st.data_cull.response=opts.response;


end

function   total=noise_fit_fun (b,x,fitobj) 
%cof_names={'freq_n','amp_shot_n','back_n'};    
freq_noise_comp=b(1)*abs(lorentzian_function_1d(x,fitobj.Coefficients{'gamma','Estimate'},...
                                fitobj.Coefficients{'mu','Estimate'},...
                                fitobj.Coefficients{'amp','Estimate'},'norm','amp','derivative',1));
int_noise_comp=abs(b(2))*sqrt(voigt_function_1d(x,fitobj.Coefficients{'sigma','Estimate'},...
                                fitobj.Coefficients{'gamma','Estimate'},fitobj.Coefficients{'mu','Estimate'},...
                                fitobj.Coefficients{'amp','Estimate'},fitobj.Coefficients{'sigma','Estimate'},...
                                'norm','amp','method','approx'));                          
total=sqrt( (freq_noise_comp).^2 +(int_noise_comp).^2 + b(3)^2 );  
                     
end


