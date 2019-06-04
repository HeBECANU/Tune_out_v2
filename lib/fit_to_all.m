function to_res=fit_to_all(anal_opts_fit_to,data)

% fit for the tune out using all the data from a run (in contrast to spliting it up scanwise)



%%

%get the mask from the output of calculate_signal 
probe_dat_mask=col_vec(data.signal.dat_mask);
to_res.fit_mask=probe_dat_mask;
to_res.num_shots=sum(probe_dat_mask);
to_res.atom_num.mean=nanmean(data.mcp_tdc.num_counts(probe_dat_mask));
to_res.atom_num.std=nanstd(data.mcp_tdc.num_counts(probe_dat_mask));
to_res.start_time=data.mcp_tdc.time_create_write(1,2);

%% get the optical frequency and the square probe trap freq
probe_freq= data.blue_probe.act.mean(probe_dat_mask');
square_trap_freq_val= data.signal.square_probe_trap_freq.val(data.signal.dat_mask);
square_trap_freq_unc=data.signal.square_probe_trap_freq.unc(data.signal.dat_mask);
%

%% define the color for each shot on the plot
cdat=viridis(1e4);
c_cord=linspace(0,1,size(cdat,1));
shot_time=data.mcp_tdc.time_create_write(probe_dat_mask,1);
shot_time=shot_time-min(shot_time);
shot_time_scaled=shot_time/range(shot_time);
cdat=[interp1(c_cord,cdat(:,1),shot_time_scaled),...
    interp1(c_cord,cdat(:,2),shot_time_scaled),...
    interp1(c_cord,cdat(:,3),shot_time_scaled)];
          
%%
%now we do some fitting
%thresholds for CI
%sd         CI
%1          0.3174
%2          0.05
%3          2.699e-03
ci_size_disp=1-erf(anal_opts_fit_to.sigma_disp/sqrt(2));%one sd %confidence interval to display
ci_size_cut_outliers=1-erf(anal_opts_fit_to.sigma_cut_outliers/sqrt(2)); %confidence interval for cutting outliers

fprintf('Calculating Fits\n')

%set up the data input for the fit
xdat=probe_freq;
ydat=square_trap_freq_val;
yunc=square_trap_freq_unc;
wdat=1./yunc;
wdat(isnan(wdat))=1e-20; %if nan then set to smallest value you can
wdat=wdat./sum(wdat);
probe_freq_offset=nanmean(xdat); 

%freq_offset=362868.029*1e9; manual set for the freq_offset;

to_res.freq_offset=probe_freq_offset;
xdat=xdat-probe_freq_offset; %fits work better when data is scaled away from large offsets
xdat=xdat*anal_opts_fit_to.scale_x;

% why is this a cell seems silly
xywdat=num2cell([xdat,ydat,wdat],2);

for ii=1:2 %iterate over linear and quadratic fits
    % for a given fit order we fit twice first to determine outliers and then again without outliers to determine the TO
    
    %% first fit to all data
    fit_out=fit_poly_with_int(xywdat,ii,1,0);
    to_res.fit_all.model{ii}=fit_out.fit_mdl;
    to_res.fit_all.scale_x{ii}=anal_opts_fit_to.scale_x;
    to_res.fit_all.to_freq(ii).val=fit_out.x_intercept.val/anal_opts_fit_to.scale_x+probe_freq_offset; %give the result back in hz
    to_res.fit_all.to_freq(ii).unc=fit_out.x_intercept.unc.with_cov/anal_opts_fit_to.scale_x;%give the unc back in hz

    xsamp=col_vec(linspace(min(xdat),max(xdat),1e3)); %sample for the model curve
    [ysamp,yci]=predict(fit_out.fit_mdl,xsamp,'Prediction','observation','Alpha',ci_size_disp); %note the observation CI
    %now make a new prediction with the model but with the CI to cut out outliers
    [~,yci_cull_lim]=predict(fit_out.fit_mdl,xdat,'Prediction','observation','Alpha',ci_size_cut_outliers);
    is_outlier_mask=ydat>yci_cull_lim(:,1) & ydat<yci_cull_lim(:,2);
    
    %% plot the fit and label the outliers
    %now plot the data and the model together
    if ii==1
        plot_name='TO fits lin';
        plot_title='Tune-out fit (Linear)';
    else
        plot_name='TO fits quad';
        plot_title='Tune-out fit (Quadratic)';
     end
    
    
    stfig(plot_name,'add_stack',1);
    subplot(3,1,1)
    plot(xsamp,ysamp,'k-')
    hold on
    plot(xsamp,yci,'r-')
    scatter(xdat,ydat,30,cdat,'square','filled')
    colormap(viridis(1000))
    c =colorbar;
    c.Label.String = 'time (H)';
    caxis([0,range(shot_time)/(60*60)])
    %plot(xdat,ydat,'bx')
    xlabel(sprintf('probe beam set freq - %.3f(GHz)',probe_freq_offset*1e-9))
    ylabel('Response (Hz^2)')
    title('Good Data')
    first_plot_lims=[get(gca,'xlim');get(gca,'ylim')];
    %color the ones that will be removed
    plot(xdat(~is_outlier_mask),ydat(~is_outlier_mask),'r.','markersize',15)
    hold off

    %% repeat the fit with the non outlier data
    culed_xywdat=xywdat(is_outlier_mask);
    cdat_culled=cdat(is_outlier_mask,:);
    xydatmat=cell2mat(culed_xywdat);
    xdat_culled= xydatmat(:,1);
    ydat_culled= xydatmat(:,2);
    wdat_culled= xydatmat(:,3);
    yunc_culled=yunc(is_outlier_mask);
    shot_time_culled=shot_time(is_outlier_mask);
    
    fit_out=fit_poly_with_int(culed_xywdat,ii,0,1);
    to_res.fit_trimmed.model{ii}=fit_out.fit_mdl;
    to_res.fit_trimmed.scale_x{ii}=anal_opts_fit_to.scale_x;
    to_res.fit_trimmed.to_freq(ii).val=fit_out.x_intercept.val/anal_opts_fit_to.scale_x+probe_freq_offset; %give the result back in hz
    to_res.fit_trimmed.to_freq(ii).unc=fit_out.x_intercept.unc.with_cov/anal_opts_fit_to.scale_x;%give the unc back in hz
    mdl_culled=fit_out.fit_mdl;
    to_res.fit_trimmed.slope.val(ii)=mdl_culled.Coefficients.Estimate(2)*anal_opts_fit_to.scale_x;
    to_res.fit_trimmed.slope.unc(ii)=mdl_culled.Coefficients.SE(2)*anal_opts_fit_to.scale_x;
    

    mdl_coef=fit_out.fit_mdl.Coefficients;
    covm=fit_out.fit_mdl.CoefficientCovariance;
    if abs((covm(1,2)-covm(2,1))/mean([covm(1,2),covm(2,1)]))>1e-9
        fprintf('error covariance matrix non symeteric')
    end
    
    %% plot the second fit
    xsamp_culled=linspace(min(xdat_culled),max(xdat_culled),1e3)';
    [ysamp_culled,yci_culled]=predict(to_res.fit_trimmed.model{ii},xsamp_culled,'Alpha',ci_size_disp); %'Prediction','observation'
    %now plot the remaining data along with the fit model and the model CI
    
    title(plot_title)
    subplot(3,1,2)
    plot(xsamp_culled,ysamp_culled,'k-')
    hold on
    plot(xsamp_culled,yci_culled,'r-')
    scatter(xdat_culled,ydat_culled,30,cdat_culled,'square','filled')
    %colormap(viridis(1000))
    c =colorbar;
    c.Label.String = 'time (H)';
    caxis([min(shot_time),max(shot_time)]/(60*60))
    hold off
    xlabel(sprintf('probe beam set freq - %.3f (GHz)',probe_freq_offset*1e-9))
    %ylabel('Response (Hz^2)')
    ylabel('\omega_{trap+probe}^2-\omega_{trap}^2 (Hz^2)')
    title('Fit Outliers Removed')
    set(gca,'xlim',first_plot_lims(1,:))
    set(gca,'ylim',first_plot_lims(2,:))
    
    yfit_dat=predict(to_res.fit_trimmed.model{ii},xdat_culled,'Alpha',ci_size_disp); %'Prediction','observation'
    subplot(3,1,3)
        errorbar(xdat_culled,ydat_culled-yfit_dat,yunc_culled,'o','CapSize',4,'MarkerSize',1,'Color',[1,1,1]*0.5,...
    'MarkerFaceColor',[1,1,1]*0.6,'LineWidth',1.5)
    hold on
    scatter(xdat_culled,ydat_culled-yfit_dat,30,cdat_culled,'square','filled')
    hold off 
    xlabel(sprintf('probe beam set freq - %.3f (GHz)',probe_freq_offset*1e-9))
    %ylabel('Response (Hz^2)')
    ylabel('Residduals (Hz^2)')
    
    
    set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1600, 900])
    saveas(gcf,[anal_opts_fit_to.global.out_dir,plot_name,'.png'])
    saveas(gcf,[anal_opts_fit_to.global.out_dir,plot_name,'.fig'])
    
    
    
    %% estimate the error in the x intercept with a bootstrap of the fit
    % use a wrapper function to return just the x intercept value
    boot_fun=@(x) wrap_fit_for_boot(x,ii,anal_opts_fit_to.scale_x);
    boot=bootstrap_se(boot_fun,culed_xywdat,...
        'plots',true,...
        'replace',true,...
        'samp_frac_lims',[0.05,1],...%[0.005,0.9]
        'num_samp_frac',10,... %20
        'num_samp_rep',3e1,... %1e2
        'plot_fig_name','TO fit bootstrap',...
        'save_multi_out',0);

     to_res.fit_trimmed.to_freq(ii).unc
     
    to_res.fit_trimmed.boot{ii}=boot;
    to_res.fit_trimmed.to_unc_boot{ii}=boot.results.se_fun_whole;
    to_res.fit_trimmed.single_shot_unc{ii}=to_res.fit_trimmed.to_unc_boot{ii}*sqrt(numel(xdat_culled));

    %find the single shot confidence interval (might need to edit this a
    %bit)
    if ii==1
        cross_xval=-mdl_culled.Coefficients.Estimate(1)/mdl_culled.Coefficients.Estimate(2);
        [cross_yval,cross_yci]=predict(mdl_culled,cross_xval,'Prediction','observation','Alpha',1-erf(1/sqrt(2)));
        if abs(cross_yval)>1e-2, error('not crossing zero here') ,end
        cross_yci=diff(cross_yci)/2;
        to_res.fit_trimmed.single_shot_uncert_fit{ii}=abs(cross_yci*(1/mdl_culled.Coefficients.Estimate(2)))/anal_opts_fit_to.scale_x;
    end
        
    %normalize by the CI at the TO
    stfig(sprintf('Senistivity Graph (order=%u)',ii),'add_stack',1);
    clf
    set(gcf,'color','w')
    plot(xsamp_culled,ysamp_culled/cross_yci,'k-')
    hold on
    plot(xsamp_culled,yci_culled/cross_yci,'r-')
    plot(xdat_culled,ydat_culled/cross_yci,'bx')
    hold off
    xlabel(sprintf('probe beam set freq - %.3f (GHz)',probe_freq_offset*1e-9))
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
temp_to_val=cell2mat({to_res.fit_trimmed.to_freq(:).val});
stfig('order dependence','add_stack',1);
errorbar(1:2,temp_to_val./1e6-mean(temp_to_val)*1e-6,cell2mat(to_res.fit_trimmed.to_unc_boot)./1e6) %with culled data
set(gcf,'color','w')
xlabel('Order of fit')
ylabel('Tune-out value (MHz)')
ylabel(sprintf('Tune-out value - %.3f (MHz)',mean(temp_to_val)*1e-6))
title('Sensitivity to fit Order')
xlim([0 3])
%% Have a look at the residuals
%res = ydat_culled-predict(mdl_culled,xdat_culled','Alpha',0.2)';


%% Make a clean plot with single markers per freq

plot_offset=(to_res.fit_trimmed.to_freq(1).val-probe_freq_offset)*1e-9;

bin_center=uniquetol(xdat_culled,1e-3);
bin_edges=cat(1,-inf,(bin_center(2:end)+bin_center(1:end-1))/2,inf);
num_bin=numel(bin_edges)-1;
group_x_mean = ones(1,num_bin);
group_y_stat=nan(numel(num_bin),3);
x_lims = ones(2,num_bin);
for jj=1:num_bin
    mask=xdat_culled>bin_edges(jj) & xdat_culled<bin_edges(jj+1);
    group_x_mean(jj) = nanmean(xdat_culled(mask));
    group_y_stat(jj,1)=sum(ydat_culled(mask).*wdat_culled(mask))./sum(wdat_culled(mask));
    group_y_stat(jj,2)=sqrt(nanvar(ydat_culled(mask),wdat_culled(mask)));
    group_y_stat(jj,3)=sewm(ydat_culled(mask),wdat_culled(mask)/sum(wdat_culled(mask)));
    try
        x_lims(:,jj) = [min(xdat_culled(mask)); max(xdat_culled(mask))];
    catch
        x_lims(:,jj) = [0;0];
    end
end

y_sd = (group_y_stat(:,2));
y_se = (group_y_stat(:,3));
xneg = group_x_mean - x_lims(1,:);
xpos = x_lims(2,:) - group_x_mean;

%Finally plot a nice version of the quad fit
%set up the colors to use
colors_main=[[233,87,0];[33,188,44];[0,165,166]];
%colors_main = [[75,151,201];[193,114,66];[87,157,95]];
font_name='cmr10';
font_size_global=14;
%we add another offset so that the plot is about the TO
group_x_mean=group_x_mean-plot_offset;

colors_main=colors_main./255;
lch=colorspace('RGB->LCH',colors_main(:,:));
lch(:,1)=lch(:,1)+20;
colors_detail=colorspace('LCH->RGB',lch);
%would prefer to use srgb_2_Jab here
color_shaded=colorspace('RGB->LCH',colors_main(3,:));
color_shaded(1)=125;
color_shaded=colorspace('LCH->RGB',color_shaded);

stfig('nice to fig','add_stack',1);
clf
[ysamp_culled,yci_culled]=predict(to_res.fit_trimmed.model{1},xsamp_culled,'Prediction','observation','Alpha',ci_size_disp); %'Prediction','observation'
patch([xsamp_culled', fliplr(xsamp_culled')]-plot_offset, [yci_culled(:,1)', fliplr(yci_culled(:,2)')], color_shaded,'EdgeColor','none');  %
hold on
xl=xlim;
line(xl,[0,0],'color','k','LineWidth',1)
plot(xsamp_culled-plot_offset,yci_culled','r','color',colors_main(3,:),'LineWidth',1.5);
plot(xsamp_culled-plot_offset,ysamp_culled,'-','color',colors_main(2,:),'LineWidth',1.5)
errorbar(group_x_mean,group_y_stat(:,1),y_sd...
     ,'o','CapSize',0,'Marker','none','Color',colors_detail(1,:),...
     'LineWidth',1.5);
errorbar(group_x_mean,group_y_stat(:,1),[],[],...
      xneg,xpos,'o','CapSize',0,'Marker','none','Color',colors_detail(1,:),...
      'LineWidth',1.5);
errorbar(group_x_mean,group_y_stat(:,1),y_se,...
    'o','CapSize',0,'MarkerSize',5,'Color',colors_main(1,:),...
    'MarkerFaceColor',colors_detail(1,:),'LineWidth',2.5);
% 
% errorbar(group_x_mean,ydat_chunks(:,1),ydat_chunks(:,2),'o','CapSize',0,'MarkerSize',5,'Color',colors_main(1,:),...
%     'MarkerFaceColor',colors_detail(1,:),'LineWidth',1.5);
hold off
xlim(xl)
xlabel(sprintf('Probe Beam (Optical) Frequency - %.3f (GHz)',probe_freq_offset*1e-9-plot_offset)) %TODO the pobe beam should be converted to the blue size much earlier in code 
ylabel('Response (Hz^2)')
set(gca,'xlim',[floor(min(group_x_mean)),ceil(max(group_x_mean))])
set(gca,'ylim',first_plot_lims(2,:))
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1600, 900])
set(gcf,'color','w')
set(gca,'FontSize',font_size_global,'FontName',font_name)
plot_name='nice_to_fig';
saveas(gcf,[anal_opts_fit_to.global.out_dir,plot_name,'.png'])
saveas(gcf,[anal_opts_fit_to.global.out_dir,plot_name,'.fig'])

end


function out=wrap_fit_for_boot(xywdat,ii,scale)
    fit_out=fit_poly_with_int(xywdat,ii,0,1);
    out=fit_out.x_intercept.val/scale;
end


