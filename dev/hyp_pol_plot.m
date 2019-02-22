%Script that scrapes the analysed data from dirs (currently messy but works)
clear all
%% setup directories
loop_config.dir = {'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20190213_to_hwp_168.5_nuller_reconfig_pdset_0.4v\';
    'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20190213_to_hwp_168.5_nuller_reconfig_pdset_0.7v\';
    'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20190213_to_hwp_168.5_nuller_reconfig_pdset_0.8v\';
    'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20190214_to_hwp_168.5_nuller_reconfig_new_fiber_pdset_1.0v\';
    'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20190214_to_hwp_168.5_nuller_reconfig_new_fiber_pdset_2.0v\';
    'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20190214_to_hwp_168.5_nuller_reconfig_new_fiber_pdset_2.5v\';
    'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20190214_to_hwp_168.5_nuller_reconfig_pdset_0.6v\'
    };

%% Get data
data = to_data(loop_config);

% Extract variables to plot
grads = data.drift.grad{:,1}./2;
grads_unc = data.drift.grad{:,2}./2;
to_vals = data.drift.to_val{:,1};
to_vals_unc = data.drift.to_val{:,2};

grads_avg = data.main.grad./2;
to_vals_run = data.main.lin_fit{:,1};
to_vals_run_unc =  data.main.lin_fit{:,2};

disp_config.num_bin = 12; 
disp_config.colors_main = [[233,87,0];[33,188,44];[0,165,166]];
disp_config.plot_title = 'Hyperpolarizability';
disp_config.font_name = 'cmr10';
disp_config.font_size_global=14;
disp_config.mdl_fun = @(b,x) b(1)+b(2).*x(:,1);
disp_config.beta0 = [1e14,1e5];
disp_config.opts=statset('nlinfit');
disp_config.fig_number=821;
plot_sexy(disp_config,grads,to_vals,to_vals_unc)



%%
% % Plot 'em all with errors
% sfigure(9009);
% errorbar(grads,to_vals,to_vals_unc,'kx','CapSize',0)
% xlabel('Gradient of Signal')
% ylabel(' Tune-out value (MHz)')
% set(gcf,'color','w')
% box on
% 
% sfigure(9010);
% errorbar(grads_avg,to_vals_run,to_vals_run_unc,'ko','CapSize',0)
% xlabel('Gradient of Signal')
% ylabel(' Tune-out value (MHz)')
% set(gcf,'color','w')
% box on
% %%
% % Set fit weights & normalize
% weights=1./to_vals_unc.^2';
% weights=weights/sum(weights);
% num_bin = 12;
% x_range = max(grads)-min(grads);
% bin_size = x_range/num_bin;
% c_data = viridis(num_bin);
% x_grouped = ones(1,num_bin);
% y_grouped=nan(numel(num_bin),2);
% x_lims = ones(2,num_bin);
% for jj=0:(num_bin-1)
%     bin_centre = bin_size*0.5+jj*bin_size+min(grads);
%     x_mask = (abs(grads-bin_centre)<(bin_size*0.5));
%     x_grouped(jj+1) = nanmean(grads(x_mask));
%     y_grouped(jj+1,1)=sum(to_vals(x_mask).*weights(x_mask)')./sum(weights(x_mask));
%     y_grouped(jj+1,2)=sqrt(nanvar(to_vals(x_mask),weights(x_mask)));
%     try
%         x_lims(:,jj+1) = [min(grads(x_mask)); max(grads(x_mask))];
%     catch
%         x_lims(:,jj+1) = [0;0];
%     end
% end
% 
% yneg = (y_grouped(:,2).*1e-6)./2;
% ypos = (y_grouped(:,2).*1e-6)./2;
% xneg = x_grouped - x_lims(1,:);
% xpos = x_lims(2,:) - x_grouped;
% 
% %Finally plot a nice version of the quad fit
% %set up the colors to use
% colors_main=[[233,87,0];[33,188,44];[0,165,166]];
% plot_title='Hyperpolarizability';
% font_name='cmr10';
% font_size_global=14;
% 
% 
% colors_main=colors_main./255;
% lch=colorspace('RGB->LCH',colors_main(:,:));
% lch(:,1)=lch(:,1)+20;
% colors_detail=colorspace('LCH->RGB',lch);
% %would prefer to use srgb_2_Jab here
% color_shaded=colorspace('RGB->LCH',colors_main(3,:));
% color_shaded(1)=125;
% color_shaded=colorspace('LCH->RGB',color_shaded);
% 
% mdl_fun = @(b,x) b(1)+b(2).*x(:,1);
% beta0 = [1e14,1e5];
% opts = statset('nlinfit');
% 
% fit_mdl = fitnlm(grads',to_vals',mdl_fun,beta0,'Options',opts,'Weight',weights);%'ErrorModel','combined'
% ci_size_disp = 1-erf(1/sqrt(2));
% %%
% sfigure(822);
% clf
% x_grouped_pad = linspace(-1, x_grouped(end)*2.01,6);
% [ysamp_culled,yci_culled]=predict(fit_mdl,x_grouped_pad','Alpha',ci_size_disp); %'Prediction','observation'
% %we add another offset so that the plot is about the TO
% plot_offset=predict(fit_mdl,0);
% patch([x_grouped_pad, fliplr(x_grouped_pad)], ([yci_culled(:,1)', fliplr(yci_culled(:,2)')]-plot_offset).*1e-9, color_shaded,'EdgeColor','none');  %
% hold on
% plot(x_grouped_pad,(yci_culled'-plot_offset).*1e-6,'r','color',colors_main(3,:),'LineWidth',1.5);
% xl=xlim;
% %line(xl,[0,0],'color','k','LineWidth',1)
% title(plot_title)
% plot(x_grouped_pad,(ysamp_culled-plot_offset).*1e-6,'-','color',colors_main(2,:),'LineWidth',1.5)
% 
% errorbar(x_grouped,(y_grouped(:,1)-plot_offset).*1e-6,yneg,ypos,xneg,xpos,'o','CapSize',0,'MarkerSize',5,'Color',colors_main(1,:),...
%     'MarkerFaceColor',colors_detail(1,:),'LineWidth',1.5);
% xlabel('Gradient of Signal (Hz^2/Hz) \times 10^{-9}')
% ylabel(sprintf('Tune-out value - %.3f (MHz)',plot_offset.*1e-6))
% set(gca,'xlim',[0,ceil(max(x_grouped)*1.05)])
% %set(gca,'ylim',first_plot_lims(2,:))
% set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1600, 900])
% set(gcf,'color','w')
% set(gca,'FontSize',font_size_global,'FontName',font_name)
% box on