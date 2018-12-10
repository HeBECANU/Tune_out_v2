function to_res=fit_to_nonlin(anal_opts_fit_to,data)

temp_cal=data.mcp_tdc.probe.calibration'; %because its used a lot make a temp var for calibration logic vector
temp_cal(isnan(temp_cal))=1;    
%manual bootstrap rand(size(data.osc_fit.ok.rmse))>0.9
probe_dat_mask=data.osc_fit.ok.all & ~temp_cal &  ~isnan(data.wm_log.proc.probe.freq.act.mean)'...
    & ~isnan(data.osc_fit.trap_freq_recons);

to_res.num_shots=sum(probe_dat_mask);
to_res.fit_mask=probe_dat_mask;

probe_freq= data.wm_log.proc.probe.freq.act.mean(probe_dat_mask')*1e6;
trap_freq=data.osc_fit.trap_freq_recons(probe_dat_mask)';
trap_freq_unc=data.osc_fit.trap_freq_recons_unc(probe_dat_mask)';
cal_trap_freq=data.cal.freq_drift_model(data.mcp_tdc.time_create_write(probe_dat_mask,1));
cal_trap_freq_unc=data.cal.unc;
delta_trap_freq=trap_freq-cal_trap_freq;
square_trap_freq= (trap_freq).^2-(cal_trap_freq).^2;
square_trap_freq_unc=sqrt(2).*sqrt((trap_freq_unc.*trap_freq).^2+(cal_trap_freq_unc.*cal_trap_freq).^2);
%define the color for each shot on the plot
cdat=viridis(1e4);
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
    ylabel('square difference in freq (\omega_{Net}^2-\omega_{cal}^2)')
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
ci_size_disp=1-erf(1/sqrt(2));%one sd %confidence interval to display
ci_size_cut_outliers=1-erf(2.5/sqrt(2));%1-erf(zvalue/sqrt(2)) %confidence interval for cutting outliers

fprintf('Calculating Fits\n')
%select the data in some freq range and that has an ok number

%set up the data input for the fit
xdat=probe_freq(~isnan(probe_freq))';
ydat=square_trap_freq(~isnan(probe_freq))';
wdat=1./square_trap_freq_unc(~isnan(probe_freq))';
wdat(isnan(wdat))=1e-20; %if nan then set to smallest value you can
wdat=wdat./sum(wdat);
cdat=cdat(~isnan(probe_freq),:);
if exist('new_to_freq_val','var') %if the TO has been calculated before use that as the center
    freq_offset=new_to_freq_val;
else
    freq_offset=nanmean(xdat); %otherwise use the center of the range
end
%freq_offset=362868.029*1e9; manual set for the freq_offset;
to_res.freq_offset=freq_offset;
xdat=xdat-freq_offset; %fits work better when thery are scaled reasonably
xdat=xdat*anal_opts_fit_to.scale_x;
xywdat=num2cell([xdat;ydat;wdat],1);

for ii=1:2 %iterate over linear and quadratic fits

    [to_freq,mdl_all]=fit_poly_data(xywdat,ii);
    to_res.fit_all.model{ii}=mdl_all;

    to_res.fit_all.to_freq{ii}=to_freq/anal_opts_fit_to.scale_x+freq_offset;

    xsamp=linspace(min(xdat),max(xdat),1e3)'; %sample for the model curve
    [ysamp,yci]=predict(mdl_all,xsamp,'Prediction','observation','Alpha',ci_size_disp); %note the observation CI
    %now make a new prediction with the model but with the CI to cut out outliers
    [~,yci_cull_lim]=predict(mdl_all,xdat','Prediction','observation','Alpha',ci_size_cut_outliers);
    is_outlier_idx=ydat>yci_cull_lim(:,1)' & ydat<yci_cull_lim(:,2)';

    %now plot the data and the model together
    sfigure(661+ii);
    set(gcf,'color','w')
    subplot(1,2,1)
    plot(xsamp,ysamp,'k-')
    hold on
    plot(xsamp,yci,'r-')
    scatter(xdat,ydat,30,cdat,'square','filled')
    colormap(viridis(1000))
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
    wdat_culled=wdat(is_outlier_idx);
    cdat_culled=cdat(is_outlier_idx,:);

    culed_xywdat=num2cell([xdat_culled;ydat_culled;wdat_culled],1);

    [to_freq,mdl_culled]=fit_poly_data(culed_xywdat,ii);

    to_res.fit_trimmed.to_freq{ii}=to_freq/anal_opts_fit_to.scale_x+freq_offset;
    to_res.fit_trimmed.model{ii}=mdl_culled;
    mdl_coef=mdl_culled.Coefficients;
    covm=mdl_culled.CoefficientCovariance;
    if abs((covm(1,2)-covm(2,1))/mean([covm(1,2),covm(2,1)]))>1e-9
        fprintf('error covariance matrix non symeteric')
    end
    if ii==1 %propogation for the linear intercept
        to_res.fit_trimmed.to_unc_fit_no_cov{ii}=abs((mdl_coef.Estimate(1)/mdl_coef.Estimate(2)))*...
            sqrt(...
            (mdl_coef.SE(1)/mdl_coef.Estimate(1))^2+...
            (mdl_coef.SE(2)/mdl_coef.Estimate(2))^2 ...
            )/anal_opts_fit_to.scale_x;

        %calculate the error propagation with the covariance included
        %generaly makes a tiny change ~1e-3
        to_res.fit_trimmed.to_unc_fit{ii}=abs((mdl_coef.Estimate(1)/mdl_coef.Estimate(2)))*...
            sqrt(...
            (mdl_coef.SE(1)/mdl_coef.Estimate(1))^2+...
            (mdl_coef.SE(2)/mdl_coef.Estimate(2))^2 ...
            - 2*covm(2,1)/(mdl_coef.Estimate(1)*mdl_coef.Estimate(2))...
            )/anal_opts_fit_to.scale_x;
    else %propogation for quadratic intercept
        det = sqrt(mdl_coef.Estimate(2)^2-4*mdl_coef.Estimate(1)*mdl_coef.Estimate(3));
        c = mdl_coef.Estimate(1);
        b = mdl_coef.Estimate(2);
        a= mdl_coef.Estimate(3);
        unc_a_c_2 = (mdl_coef.SE(1)^2+(2*c*a+b*(det-b))^2/(4*a^4)*mdl_coef.SE(3)^2)*1/det^2;
        unc_b_2 = (b/det-1)^2*1/(4*a^2)*mdl_coef.SE(2)^2;
        to_res.fit_trimmed.to_unc_fit_no_cov{ii}=sqrt(...
            unc_a_c_2+...
            unc_b_2 )/anal_opts_fit_to.scale_x;

        %calculate the error propagation with the covariance included
        %generaly makes a tiny change ~1e-3
        unc_ab = (-b^3+3*a*b*c+det*(b^2-a*c))*covm(3,2);
        unc_ac = (a*b^2-2*a^2*c-a*b*det)*covm(3,1);
        unc_bc = (-a^2*b+a^2*det)*covm(2,1);
        to_res.fit_trimmed.to_unc_fit{ii}=sqrt(to_res.fit_trimmed.to_unc_fit_no_cov{ii}^2+...
        1/(a^3*det^2)*(unc_ab+unc_ac+unc_bc)...
            /anal_opts_fit_to.scale_x^2);
    end
    poly_mdl = @(in) fit_poly_data(in,ii); 
    boot=bootstrap_se(poly_mdl,culed_xywdat,...
        'plots',true,...
        'replace',true,...
        'samp_frac_lims',[0.05,1],...%[0.005,0.9]
        'num_samp_frac',10,...
        'num_samp_rep',3e1,...
        'plot_fig_num',89,...
        'save_multi_out',0);

    to_res.fit_trimmed.boot{ii}=boot;
    to_res.fit_trimmed.to_unc_boot{ii}=boot.se_opp/anal_opts_fit_to.scale_x;
    to_res.fit_trimmed.single_shot_uncert_boot{ii}=to_res.fit_trimmed.to_unc_boot{ii}*sqrt(numel(xdat_culled));

    xsamp_culled=linspace(min(xdat_culled),max(xdat_culled),1e3)';
    [ysamp_culled,yci_culled]=predict(mdl_culled,xsamp_culled,'Alpha',0.2); %'Prediction','observation'
    %now plot the remaining data along with the fit model and the model CI
    if ii==1
        plot_name='TO_fits_lin';
        plot_title='Tune-out fit (Linear)';
    else
        plot_name='TO_fits_quad';
        plot_title='Tune-out fit (Quadratic)';
    end
    sfigure(661+ii);
    title(plot_title)
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
    %ylabel('Response (Hz^2)')
    ylabel('\omega_{trap+probe}^2-\omega_{trap}^2 (Hz^2)')
    title('Fit Outliers Removed')
    set(gca,'xlim',first_plot_lims(1,:))
    set(gca,'ylim',first_plot_lims(2,:))
    
    set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1600, 900])
    saveas(gcf,[anal_opts_fit_to.global.out_dir,plot_name,'.png'])
    saveas(gcf,[anal_opts_fit_to.global.out_dir,plot_name,'.fig'])

    %find the single shot confidence interval (might need to edit this a
    %bit)
    if ii==1
        cross_xval=-mdl_culled.Coefficients.Estimate(1)/mdl_culled.Coefficients.Estimate(2);
        [cross_yval,cross_yci]=predict(mdl_culled,cross_xval,'Prediction','observation','Alpha',0.3174);
        if abs(cross_yval)>1e-2, error('not crossing zero here') ,end
        cross_yci=diff(cross_yci)/2;
        to_res.fit_trimmed.single_shot_uncert_fit{ii}=abs(cross_yci*(1/mdl_culled.Coefficients.Estimate(2)))/anal_opts_fit_to.scale_x;
    end
        
    %normalize by the CI at the TO
    figure(761+ii);
    clf
    set(gcf,'color','w')
    plot(xsamp_culled,ysamp_culled/cross_yci,'k-')
    hold on
    plot(xsamp_culled,yci_culled/cross_yci,'r-')
    plot(xdat_culled,ydat_culled/cross_yci,'bx')
    hold off
    xlabel(sprintf('probe beam set freq - %.3f (GHz)',freq_offset*1e-9))
    ylabel('Response scaled to sample SD')
    if ii==1
        title_str='Senistivity Graph (Linear)';
        plot_name='Sens_graph_lin';
    else
        title_str='Senistivity Graph (Quadratic)';
        plot_name='Sens_graph_quad';
    end
    title(title_str)
    set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1600, 900])
    saveas(gcf,[anal_opts_fit_to.global.out_dir,plot_name,'.png'])
    saveas(gcf,[anal_opts_fit_to.global.out_dir,plot_name,'.fig'])
end
%% See how the intercept error changes with order of fit (compare TO values)
%scatter(1:4,cell2mat(to_res.fit_all.to_freq)) %to with all the data
temp_to=cell2mat(to_res.fit_trimmed.to_freq);
sfigure(1234);
errorbar(1:2,temp_to./1e6-mean(temp_to)*1e-6,cell2mat(to_res.fit_trimmed.to_unc_boot)./1e6) %with culled data
set(gcf,'color','w')
xlabel('Order of fit')
ylabel('Tune-out value (MHz)')
ylabel(sprintf('Tune-out value - %.3f (MHz)',mean(temp_to)*1e-6))
title('Sensitivity to fit Order')
xlim([0 3])
%% Have a look at the residuals
res = ydat_culled-predict(mdl_culled,xdat_culled','Alpha',0.2)';
num_bin = 20;
x_range = max(xdat_culled)-min(xdat_culled);
bin_size = x_range/num_bin;
c_data = viridis(num_bin);
x_res_grouped = ones(1,num_bin);
y_res_grouped = cell(1,num_bin);
ydat_chunks=nan(numel(num_bin),2);
sfigure(487);
clf
for jj=0:(num_bin-1)
%     hold on
    bin_centre = bin_size*0.5+jj*bin_size+min(xdat_culled);
    x_res_grouped(jj+1) = bin_centre;
    x_mask = (abs(xdat_culled-bin_centre)<(bin_size*0.5));
    y_res_grouped{:,jj+1} = res(x_mask)';
    ydat_chunks(jj+1,1)=nanmean(ydat_culled(x_mask));
    ydat_chunks(jj+1,2)=nanstd(ydat_culled(x_mask));
end
violin(y_res_grouped,'x',x_res_grouped,'facecolor',c_data,'edgecolor','none','bw',10,'mc','k','medc','r-.');
ylabel('Residuals in Signal')
xlabel(sprintf('Tune-out value - %.3f (MHz)',freq_offset*1e-6))
set(gcf,'color','w')
sfigure(476);
histogram(res,'FaceAlpha',0.45)
xlabel('Residuals in Signal')
ylabel('Count')
set(gcf,'color','w')


%Finally plot a nice version of the quad fit
%set up the colors to use
colors_main=[[233,87,0];[33,188,44];[0,165,166]];
plot_title='Tune-Out Fit';
font_name='cmr10';
font_size_global=14;

%we add another offset so that the plot is about the TO
plot_offset=1*(to_res.fit_trimmed.to_freq{2}-freq_offset)*1e-9;
x_res_grouped=x_res_grouped-plot_offset;


colors_main=colors_main./255;
lch=colorspace('RGB->LCH',colors_main(:,:));
lch(:,1)=lch(:,1)+20;
colors_detail=colorspace('LCH->RGB',lch);
%would prefer to use srgb_2_Jab here
color_shaded=colorspace('RGB->LCH',colors_main(3,:));
color_shaded(1)=125;
color_shaded=colorspace('LCH->RGB',color_shaded);

sfigure(82);
clf
[ysamp_culled,yci_culled]=predict(mdl_culled,xsamp_culled,'Prediction','observation','Alpha',ci_size_disp); %'Prediction','observation'
patch([xsamp_culled', fliplr(xsamp_culled')]-plot_offset, [yci_culled(:,1)', fliplr(yci_culled(:,2)')], color_shaded,'EdgeColor','none');  %
hold on
plot(xsamp_culled-plot_offset,yci_culled','r','color',colors_main(3,:),'LineWidth',1.5);
xl=xlim;
line(xl,[0,0],'color','k','LineWidth',1)
title(plot_title)
plot(xsamp_culled-plot_offset,ysamp_culled,'-','color',colors_main(2,:),'LineWidth',1.5)
errorbar(x_res_grouped,ydat_chunks(:,1),ydat_chunks(:,2),'o','CapSize',0,'MarkerSize',5,'Color',colors_main(1,:),...
    'MarkerFaceColor',colors_detail(1,:),'LineWidth',1.5);
xlabel(sprintf('Probe Beam (Optical) Frequency - %.3f (GHz)',freq_offset*1e-9-plot_offset))
ylabel('Response (Hz^2)')
set(gca,'xlim',[floor(min(x_res_grouped)),ceil(max(x_res_grouped))])
set(gca,'ylim',first_plot_lims(2,:))
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1600, 900])
set(gcf,'color','w')
set(gca,'FontSize',font_size_global,'FontName',font_name)
end

function [to_freq,fit_mdl]=fit_poly_data(in,n)
inmat=cell2mat(in);
xdat=inmat(1,:);
ydat=inmat(2,:);
wdat=inmat(3,:);

modelfun = @(b,x) (repmat(x(:,1),1,n+1).^(0:n))*(b(1:(n+1))'); %simple linear model
opts = statset('nlinfit');
%opts.RobustWgtFun = 'welsch' ; %a bit of robust fitting
%opts.Tune = 1;
beta0 = [1e-5,1e-2.*ones(1,n)]; %intial guesses
fit_mdl = fitnlm(xdat,ydat,modelfun,beta0,'Options',opts,'ErrorModel','combined'); %constant, 'Weight',wdat add for weighted fits, haven't got them quite working yet

%fzero(@(x) predict(fit_mdl,x),0)
mdl_zeros = roots(fliplr(fit_mdl.Coefficients.Estimate(:)'));
to_freq_estimate=-fit_mdl.Coefficients.Estimate(1)/fit_mdl.Coefficients.Estimate(2);
[c, indx] = min(abs(mdl_zeros-to_freq_estimate));
to_freq = mdl_zeros(indx);
end
