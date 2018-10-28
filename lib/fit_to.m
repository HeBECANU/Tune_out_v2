function to_res=fit_to(anal_opts_fit_to,data)

temp_cal=data.mcp_tdc.probe.calibration';
temp_cal(isnan(temp_cal))=1;    
%manual bootstrap rand(size(data.osc_fit.ok.rmse))>0.9
probe_dat_mask=data.osc_fit.ok.rmse & ~temp_cal & ~isnan(data.mcp_tdc.probe.freq.act.mean') &...
    ~isnan(data.mcp_tdc.probe.freq.act.mean)' & ~isnan(data.osc_fit.trap_freq_recons);

to_res.num_shots=sum(probe_dat_mask);
to_res.fit_mask=probe_dat_mask;

probe_freq= data.mcp_tdc.probe.freq.act.mean(probe_dat_mask')*1e6;
trap_freq=data.osc_fit.trap_freq_recons(probe_dat_mask)';
cal_trap_freq=data.cal.freq_drift_model(data.mcp_tdc.time_create_write(probe_dat_mask,1));
delta_trap_freq=trap_freq-cal_trap_freq;
square_trap_freq= (trap_freq).^2-(cal_trap_freq).^2;
%define the color for each shot on the plot
cdat=viridis(1000);
c_cord=linspace(0,1,size(cdat,1));
shot_time=data.mcp_tdc.time_create_write(probe_dat_mask,1);
shot_time=shot_time-min(shot_time);
shot_time_scaled=shot_time/range(shot_time);
cdat=[interp1(c_cord,cdat(:,1),shot_time_scaled),...
    interp1(c_cord,cdat(:,2),shot_time_scaled),...
    interp1(c_cord,cdat(:,3),shot_time_scaled)];
        
        
if anal_opts_fit_to.plot_inital
    figure(71)
    clf
    set(gcf,'color','w')
    subplot(2,1,1)

    scatter((probe_freq-nanmean(probe_freq))*1e-9,delta_trap_freq,30,cdat,'square','filled')
    colormap(viridis(1000))
    c =colorbar;
    c.Label.String = 'time (H)';
    caxis([0,range(shot_time)/(60*60)])
    xlabel('delta probe beam frequency (GHz)')
    ylabel('delta freq')
    subplot(2,1,2)

    scatter((probe_freq-nanmean(probe_freq))*1e-9,square_trap_freq,30,cdat,'square','filled')
    xlabel('delta probe beam frequency (GHz)')
    ylabel('square difference in freq')
    colormap(viridis(1000))
    c =colorbar;
    c.Label.String = 'time (H)';
    caxis([0,range(shot_time)/(60*60)])

    %set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1600, 900])
    plot_name='tuneout_time_graph'; 
    saveas(gcf,[anal_opts_fit_to.global.out_dir,plot_name,'.png'])
    saveas(gcf,[anal_opts_fit_to.global.out_dir,plot_name,'.fig'])
end


%%
%now we do some fitting
%thresholds for CI
%sd         CI
%1          0.3174
%2          0.05
%3          2.699e-03
ci_size_disp=0.3174;%one sd %confidence interval to display
ci_size_cut_outliers=0.01; %confidence interval for cutting outliers

fprintf('Calculating Fits\n')
%select the data in some freq range and that has an ok number

%set up the data input for the fit
xdat=probe_freq(~isnan(probe_freq))';
ydat=square_trap_freq(~isnan(probe_freq))';
cdat=cdat(~isnan(probe_freq));
if exist('new_to_freq_val','var') %if the TO has been calculated before use that as the center
    freq_offset=new_to_freq_val;
else
    freq_offset=nanmean(xdat); %otherwise use the center of the range
end
xdat=xdat-freq_offset; %fits work better when thery are scaled reasonably
xdat=xdat*anal_opts_fit_to.scale_x;
xydat=num2cell([xdat;ydat],1);

[to_freq,mdl_all]=fit_lin_data(xydat);
to_res.fit_all.to_freq=to_freq/anal_opts_fit_to.scale_x+freq_offset;
to_res.fit_all.model=mdl_all;
to_res.freq_offset=freq_offset;

xsamp=linspace(min(xdat),max(xdat),1e3)'; %sample for the model curve
[ysamp,yci]=predict(mdl_all,xsamp,'Prediction','observation','Alpha',ci_size_disp); %note the observation CI
%now make a new prediction with the model but with the CI to cut out outliers
[~,yci_cull_lim]=predict(mdl_all,xdat','Prediction','observation','Alpha',ci_size_cut_outliers);
is_outlier_idx=ydat>yci_cull_lim(:,1)' & ydat<yci_cull_lim(:,2)';

%now plot the data and the model together
sfigure(6);
set(gcf,'color','w')
subplot(1,2,1)
plot(xsamp,ysamp,'k-')
hold on
plot(xsamp,yci,'r-')
scatter(xdat,ydat,30,cdat,'square','filled')
c =colorbar;
c.Label.String = 'time (H)';
caxis([0,range(shot_time)/(60*60)])
%plot(xdat,ydat,'bx')
xlabel(sprintf('probe beam set freq - %.3f(GHz)',freq_offset*1e-9))
ylabel('Response (Hz^2)')
title('Good Data')
first_plot_lims=[get(gca,'xlim');get(gca,'ylim')];
%color the ones that will be removed
plot(xdat(~is_outlier_idx),ydat(~is_outlier_idx),'r.','markersize',15)
hold off


xdat_culled=xdat(is_outlier_idx);
ydat_culled=ydat(is_outlier_idx);
cdat_culled=cdat(is_outlier_idx);

culed_xydat=num2cell([xdat_culled;ydat_culled],1);

[to_freq,mdl_culled]=fit_lin_data(culed_xydat);
to_res.fit_trimmed.to_freq=to_freq/anal_opts_fit_to.scale_x+freq_offset;
to_res.fit_trimmed.model=mdl_culled;
mdl_coef=mdl_culled.Coefficients;
to_res.fit_trimmed.to_unc_fit=abs((mdl_coef.Estimate(1)/mdl_coef.Estimate(2)))*...
    sqrt((mdl_coef.SE(1)/mdl_coef.Estimate(1))^2+...
    (mdl_coef.SE(2)/mdl_coef.Estimate(2))^2)/anal_opts_fit_to.scale_x;


boot=bootstrap_se(@fit_lin_data,culed_xydat,...
    'plots',false,...
    'replace',true,...
    'samp_frac_lims',[0.05,0.9],...%[0.005,0.9]
    'num_samp_frac',3,...
    'num_samp_rep',1e1,...
    'plot_fig_num',89,...
    'save_multi_out',0);

to_res.fit_trimmed.boot=boot;
to_res.fit_trimmed.to_unc_boot=to_res.fit_trimmed.boot.se_opp/anal_opts_fit_to.scale_x;
to_res.fit_trimmed.single_shot_uncert_boot=to_res.fit_trimmed.to_unc_boot*sqrt(numel(xdat_culled));



xsamp_culled=linspace(min(xdat_culled),max(xdat_culled),1e3)';
[ysamp_culled,yci_culled]=predict(mdl_culled,xsamp_culled,'Alpha',0.2); %'Prediction','observation'
%now plot the remaining data along with the fit model and the model CI
sfigure(6);
subplot(1,2,2)
plot(xsamp_culled,ysamp_culled,'k-')
hold on
plot(xsamp_culled,yci_culled,'r-')
scatter(xdat_culled,ydat_culled,30,cdat_culled,'square','filled')
colormap(viridis(1000))
c =colorbar;
c.Label.String = 'time (H)';
caxis([0,range(shot_time)/(60*60)])
hold off
xlabel(sprintf('probe beam set freq - %.3f (GHz)',freq_offset*1e-9))
ylabel('Response (Hz^2)')
title('Fit Outliers Removed')
set(gca,'xlim',first_plot_lims(1,:))
set(gca,'ylim',first_plot_lims(2,:))


set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1600, 900])
plot_name='TO_fits';
saveas(gcf,[anal_opts_fit_to.global.out_dir,plot_name,'.png'])
% saveas(gcf,[anal_opts_fit_to.global.out_dir,plot_name,'.fig'])

%find the single shot confidence interval
cross_xval=-mdl_culled.Coefficients.Estimate(1)/mdl_culled.Coefficients.Estimate(2);
[cross_yval,cross_yci]=predict(mdl_culled,cross_xval,'Prediction','observation','Alpha',0.3174);
if abs(cross_yval)>1e-2, error('not crossing zero here') ,end
cross_yci=diff(cross_yci)/2;
to_res.fit_trimmed.single_shot_uncert_fit=abs(cross_yci*(1/mdl_culled.Coefficients.Estimate(2)))/anal_opts_fit_to.scale_x;

%normalize by the CI at the TO
figure(7);
clf
set(gcf,'color','w')
plot(xsamp_culled,ysamp_culled/cross_yci,'k-')
hold on
plot(xsamp_culled,yci_culled/cross_yci,'r-')
plot(xdat_culled,ydat_culled/cross_yci,'bx')
hold off
xlabel(sprintf('probe beam set freq - %.3f (GHz)',freq_offset*1e-9))
ylabel('Response scaled to sample SD')
title('Senistivity Graph ')

set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1600, 900])
plot_name='Sens_graph';
saveas(gcf,[anal_opts_fit_to.global.out_dir,plot_name,'.png'])
% saveas(gcf,[anal_opts_fit_to.global.out_dir,plot_name,'.fig'])


end

function [to_freq,fit_mdl]=fit_lin_data(in);
inmat=cell2mat(in);
xdat=inmat(1,:);
ydat=inmat(2,:);

modelfun = @(b,x) b(1)+ b(2).*x; %simple linear model
opts = statset('nlinfit');
%opts.RobustWgtFun = 'welsch' ; %a bit of robust fitting
%opts.Tune = 1;
beta0 = [1e-5,1e-2]; %intial guesses
fit_mdl = fitnlm(xdat,ydat,modelfun,beta0,'Options',opts);
to_freq=-fit_mdl.Coefficients.Estimate(1)/fit_mdl.Coefficients.Estimate(2);
end
