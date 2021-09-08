% function dum  = hyp_pol_plot()
% %Script that scrapes the analysed data from dirs (currently messy but works)
%% setup directories
loop_config.dir = {'E:\Data\Tuneout\to_main_data\20190211T1737_to_hwp_333.5_polmin_168.5_nuller_reconfig_pdset_0.2v';
    'E:\Data\Tuneout\to_main_data\20190213T1015_to_hwp_333.5_polmin_168.5_nuller_reconfig_pdset_0.8v';
    'E:\Data\Tuneout\to_main_data\20190213T1242_to_hwp_333.5_polmin_168.5_nuller_reconfig_pdset_0.7v';
    'E:\Data\Tuneout\to_main_data\20190213T2024_to_hwp_333.5_polmin_168.5_nuller_reconfig_pdset_0.4v';
    'E:\Data\Tuneout\to_main_data\20190214T1220_to_hwp_333.5_polmin_168.5_nuller_reconfig_pdset_0.6v';
    'E:\Data\Tuneout\to_main_data\20190214T1657_to_hwp_333.5_polmin_168.5_nuller_reconfig_new_fiber_pdset_2.0v';
    'E:\Data\Tuneout\to_main_data\20190214T1818_to_hwp_333.5_polmin_168.5_nuller_reconfig_new_fiber_pdset_2.5v';
    'E:\Data\Tuneout\to_main_data\20190214T1942_to_hwp_333.5_polmin_168.5_nuller_reconfig_new_fiber_pdset_1.0v'};

%% Get data
data = load_pocessed_to_data(loop_config);
%%
% Extract variables to plot
grads = 1e9*data.drift.grad.val;
grads_unc = 1e9*data.drift.grad.unc;
to_vals = data.drift.to.val.*1e-6;
to_vals_unc = data.drift.to.unc.*1e-6;


grad_cut = ~(grads<=0);

grads = grads(grad_cut);
grads_unc = grads_unc(grad_cut);
to_vals = to_vals(grad_cut);
to_vals_unc = to_vals_unc(grad_cut);

% pool them
grad_edges = 0:2:90;
grad_cens = 0.5*(grad_edges(2:end) + grad_edges(1:end-1));
num_bins = length(grad_edges)-1;
to_mean = zeros(num_bins,1);
to_se = zeros(num_bins,1);
for bidx = 1:num_bins
    mask = grads > grad_edges(bidx) & grads < grad_edges(bidx+1);
    if sum(mask)>0
        to_mean(bidx) = nanmean(to_vals(mask));
    %     to_se(bidx) = mean(to_vals_unc(mask));
        to_se(bidx) = std(to_vals(mask));
    end
end
nmask = ~isnan(to_mean);
to_mean = to_mean(nmask);
to_se = to_se(~isnan(to_se));
grad_cens = grad_cens(nmask);
% Plot 'e
% grads_avg = data.main.grad./2;
% to_vals_run = data.main.lin_fit{:,1};
% to_vals_run_unc =  data.main.lin_fit{:,2};


cli_header('- - - - - ');
cli_header('RUNNING ANALYSIS');
cli_header('- - - - - ');
% Run the various fits

alpha_cull = erf(sqrt(1/2));

% grad_pow = 2;
weights = 1./to_vals_unc.^2;
weights = weights/sum(weights);
% weights = grads.;
mdl_fun = @(b,x) b(1)+b(2).*x(:,1);
beta0 = [nanmean(to_vals),1e5];
opts = statset('MaxIter',1e4);
% full fit
cli_header('UNWT FIT ALL DATA');
[full_mdl,gof,fit_out] = fit(grads,to_vals,'poly1')
full_ci=confint(full_mdl,1-erf(1/sqrt(2)))-[full_mdl.p1,full_mdl.p2];
data_pi = predint(full_mdl,grads,1-alpha_cull,'observation','off');
outlier_z = abs(fit_out.residuals)./(0.5*diff(data_pi,[],2));
outlier_mask = outlier_z<5;


alpha = erf(sqrt(1/2));
cli_header('WEIGHTED FIT ALL DATA');
% full weighted fit
[wt_mdl,wt_gof,wt_out] = fit(grads,to_vals,'poly1','Weights',weights)%'ErrorModel','combined'
wt_ci=confint(wt_mdl,1-alpha)-[wt_mdl.p1,wt_mdl.p2];

% Masking outliers
mwts = weights(outlier_mask);
mwts = mwts/sum(mwts); % for convenience later

cli_header('WEIGHTED FIT OUTLIER CULLED');
[mask_mdl,mask_gof,mask_out] = fit(grads(outlier_mask),to_vals(outlier_mask),'poly1','Weights',mwts)
mask_mdl_full= fitlm(grads(outlier_mask),to_vals(outlier_mask),'Weights',mwts)
mask_mdl_full.ModelCriterion
% b etc = regress(x,y)
mask_ci=confint(mask_mdl,0.95)-[mask_mdl.p1,mask_mdl.p2];

cli_header('UNWT FIT OUTLIER CULLED');
% weighted culled fit
[uwt_mask_mdl,uwt_mask_gof,uwt_mask_out] = fit(grads(outlier_mask),to_vals(outlier_mask),'poly1')
uwt_mask_mdl_full = fitlm(grads(outlier_mask),to_vals(outlier_mask))
uwt_mask_ci=confint(uwt_mask_mdl,0.95)-[uwt_mask_mdl.p1,uwt_mask_mdl.p2];
uwt_mask_mdl_full.ModelCriterion

sum((mask_mdl(grads(outlier_mask))-to_vals(outlier_mask)).^2)
sum((uwt_mask_mdl(grads(outlier_mask))-to_vals(outlier_mask)).^2)


% Pool culled vals
% pool the culled vals
to_mean_cull = zeros(num_bins,1);
to_se_cull = zeros(num_bins,1);

grads_cull = grads(outlier_mask);
to_vals_cull = to_vals(outlier_mask);
to_vals_unc_cull = to_vals_unc(outlier_mask);
for bidx = 1:num_bins
    mask = grads_cull > grad_edges(bidx) & grads_cull < grad_edges(bidx+1);
    if sum(mask)>0
        to_mean_cull(bidx) = nanmean(to_vals_cull(mask));
        to_se_cull(bidx) = std(to_vals_cull(mask));
    end
end
nmask_cull = ~isnan(to_mean_cull);
to_mean_cull = to_mean_cull(nmask_cull);
to_se_cull = to_se_cull(~isnan(to_se_cull));

% % Print output
% alpha = 0.95;
cli_header('----')
cli_header('fit slope %.2f(%.2f,%.2f)',full_mdl.p1,mask_ci(:,1));
cli_header('culled fit slope %.2f(%.2f,%.2f)',uwt_mask_mdl.p1,uwt_mask_ci(:,1));
cli_header('weighted fit slope %.2f(%.2f,%.2f)',wt_mdl.p1,wt_ci(:,1));
cli_header('weighted culled fit slope %.2f(%.2f,%.2f)',mask_mdl.p1,mask_ci(:,1));
cli_header('----')
cli_header('fit intercept %.2f(%.2f,%.2f)',full_mdl.p2,mask_ci(:,2));
cli_header('culled fit intercept %.2f(%.2f,%.2f)',uwt_mask_mdl.p2,uwt_mask_ci(:,2));
cli_header('weighted fit intercept %.2f(%.2f,%.2f)',wt_mdl.p2,wt_ci(:,2));
cli_header('weighted culled fit intercept %.2f(%.2f,%.2f)',mask_mdl.p2,mask_ci(:,2));
cli_header('----')
test_val = 30;
% alpha = erf(1/sqrt(2));
% alpha = 0.95;
cli_header('Shift at %.2f Hz^2/GHz (%.2f conf):',test_val,alpha);
cli_header('Shift %.2f(%.2f,%.2f)',full_mdl(test_val)-full_mdl.p2,predint(full_mdl,test_val,alpha,'functional')-full_mdl(test_val));
cli_header('culled Shiftt %.2f(%.2f,%.2f)',uwt_mask_mdl(test_val)-uwt_mask_mdl.p2,predint(uwt_mask_mdl,test_val,alpha,'functional')-uwt_mask_mdl(test_val));
cli_header('weighted Shift %.2f(%.2f,%.2f)',wt_mdl(test_val)-wt_mdl.p2,predint(wt_mdl,test_val,alpha,'functional')-wt_mdl(test_val));
cli_header('weighted culled Shift %.2f(%.2f,%.2f)',mask_mdl(test_val)-mask_mdl.p2,predint(mask_mdl,test_val,alpha,'functional')-mask_mdl(test_val));
cli_header('----')
cli_header('Shift at %.2f Hz^2/GHz (%.2f obs):',test_val,alpha);
cli_header('Shift %.2f(%.2f,%.2f)',full_mdl(test_val)-full_mdl.p2,predint(full_mdl,test_val,alpha,'observation')-full_mdl(test_val));
cli_header('culled Shiftt %.2f(%.2f,%.2f)',uwt_mask_mdl(test_val)-uwt_mask_mdl.p2,predint(uwt_mask_mdl,test_val,alpha,'observation')-uwt_mask_mdl(test_val));
cli_header('weighted Shift %.2f(%.2f,%.2f)',wt_mdl(test_val)-wt_mdl.p2,predint(wt_mdl,test_val,alpha,'observation')-wt_mdl(test_val));
cli_header('weighted culled Shift %.2f(%.2f,%.2f)',mask_mdl(test_val)-mask_mdl.p2,predint(mask_mdl,test_val,alpha,'observation')-mask_mdl(test_val));


% Set up plot elements

plotrange = linspace(0,90);
% alpha = erf(1/sqrt(2));
% alpha = .95;
to_offset = mask_mdl.p2;
h=2;
w=1;

%set up the colors to use
colors_main=[[233,87,0];[33,188,44];[0,165,166]]/255;
lch=colorspace('RGB->LCH',colors_main(:,:));
lch(:,1)=lch(:,1)+20;
colors_detail=colorspace('LCH->RGB',lch);


% full fit
y_mdl = full_mdl(plotrange);

pred_vals=full_mdl(grads);
pred_PI = predint(full_mdl,grads,erf(1/sqrt(2)),'observation','off')-pred_vals;
pred_std = mean(abs(pred_PI),2);
chi_square = sum((fit_out.residuals./pred_std).^2);
chi_square_per_dof = sum((fit_out.residuals./pred_std).^2)/gof.dfe

% weighted fit
y_wt = wt_mdl(plotrange);
y_wt_ci = predint(wt_mdl,plotrange,1-alpha);
% Unsure how to do chi-squared for weighted fit
pred_vals=wt_mdl(grads);
pred_PI = predint(wt_mdl,grads,alpha,'observation','off')-pred_vals;
pred_std = mean(abs(pred_PI),2);
chi_square = sum((wt_out.residuals.^2.*weights));
chi_square_per_dof_weighted = chi_square/wt_gof.dfe
chi_factor = chi_square_per_dof_weighted/chi_square_per_dof;

% - seems very small, suggesting 

% outlier-culled fit
y_mask_uwt = uwt_mask_mdl(plotrange);
mdl_pi = predint(uwt_mask_mdl,plotrange,alpha,'observation','off');
mdl_ci = predint(uwt_mask_mdl,plotrange,alpha,'functional','off');
wt_mdl_ci = predint(mask_mdl,plotrange,alpha,'functional','off');
wt_mdl_pi = predint(mask_mdl,plotrange,alpha,'observation','on','Weights',mwts);
wt_mdl_syn_pi = 1 * wt_mdl_pi;


pred_PI = predint(uwt_mask_mdl,grads(outlier_mask),erf(1/sqrt(2)),'observation','off')-uwt_mask_mdl(grads(outlier_mask));
pred_std = mean(abs(pred_PI),2);
chi_square = sum((uwt_mask_out.residuals./pred_std).^2);
chi_square_per_dof_mask_unwt= sum((uwt_mask_out.residuals./pred_std).^2)/uwt_mask_gof.dfe

% weighted culled fit
y_mask = mask_mdl(plotrange);





pred_PI = predint(mask_mdl,grads(outlier_mask),erf(1/sqrt(2)),'observation','off')-mask_mdl(grads(outlier_mask));
pred_std = mean(abs(pred_PI),2);
% chi_square = sum(mwts.*mask_out.residuals.^2./pred_std.^2);
% chi_square_per_dof_mask_wt= chi_square/mask_gof.dfe
 
check_mask = grads > 25 & grads < 35;
check_grad = grads(check_mask);
to_check = to_vals(check_mask);
cli_header('Val/std/se around POI: (%.3f,%.3f,%.3f)',nanmean(to_check)-to_offset,std(to_check),nanstd(to_check)/sqrt(sum(check_mask)));

% % Make plots


col_use = 2;
stfig('Hyperpolarizability Binned');
clf
hold on



piplot=fill([plotrange,fliplr(plotrange)],[mdl_pi(:,1);flipud(mdl_pi(:,2))]'-to_offset,0.7*[1,1,1],'FaceAlpha',0.3,...
    'EdgeColor','none');

wmfit=plot(plotrange,y_wt-to_offset,'k','LineWidth',2);

all_pts=errorbar(grads(outlier_mask),to_vals(outlier_mask)-to_offset,to_vals_unc(outlier_mask),'o','MarkerSize',2,'Color',colors_main(3,:),...
    'MarkerFaceColor',colors_detail(2+1,:),'MarkerEdgeColor',colors_main(2+1,:),'LineWidth',1,...
    'CapSize',0);

errorbar(grad_cens,to_mean-to_offset,to_se,'o','MarkerSize',8,'Color',colors_main(col_use,:),...
    'MarkerFaceColor',colors_detail(col_use,:),'MarkerEdgeColor',colors_main(col_use,:),'LineWidth',4,...
    'CapSize',0)

% errorbar(grad_cens+1,to_mean_cull-to_offset,to_se_cull,'o','MarkerSize',8,'Color',colors_main(col_use+1,:),...
%     'MarkerFaceColor',colors_detail(col_use+1,:),'MarkerEdgeColor',colors_main(col_use+1,:),'LineWidth',4,...
%     'CapSize',0)
% legend([ar1,ar2,ar3,ar4],{'Unwt PI','Unwt CI','Weight PI','weight CI'})


xlabel('Gradient of Signal (Hz$^2$/GHz)')
ylabel(sprintf('$f_{TO}$ - %.f (MHz)',to_offset))
ylim([-1e4,1.2e4])
set(gca,'FontSize',20)

box on
set(gca,'LineWidth',2)

col_use = 2;
% Plot raw data

stfig('Hyperpolarizability');
clf
% subplot(h,1,1)
hold on


piplot=fill([plotrange,fliplr(plotrange)],[mdl_pi(:,1);flipud(mdl_pi(:,2))]'-to_offset,0.7*[1,1,1],'FaceAlpha',.3,...
    'EdgeColor','none');

wpiplot=fill([plotrange,fliplr(plotrange)],[wt_mdl_pi(:,1);flipud(wt_mdl_pi(:,2))]'-to_offset,0.8*[0,0,1],'FaceAlpha',.3,...
    'EdgeColor','none');
% 
% synwpiplot=fill([plotrange,fliplr(plotrange)],[wt_mdl_syn_pi(:,1);flipud(wt_mdl_syn_pi(:,2))]'-to_offset,0.8*[1,0,0],'FaceAlpha',.3,...
%     'EdgeColor','none');

wciplot=fill([plotrange,fliplr(plotrange)],[wt_mdl_ci(:,1);flipud(wt_mdl_ci(:,2))]'-to_offset,0.3*[0,.3,1],'FaceAlpha',.3,...
    'EdgeColor','none');

ciplot=fill([plotrange,fliplr(plotrange)],[mdl_ci(:,1);flipud(mdl_ci(:,2))]'-to_offset,0.3*[0,.6,1],'FaceAlpha',0.6,...
    'EdgeColor','none');


wmfit=plot(plotrange,mask_mdl(plotrange)-to_offset,'k','LineWidth',2);
% wmfit=plot(plotrange,uwt_mask_mdl(plotrange)-to_offset,'k-.','LineWidth',2);


all_pts=errorbar(grads(outlier_mask),to_vals(outlier_mask)-to_offset,to_vals_unc(outlier_mask),'o','MarkerSize',6,'Color',colors_main(col_use,:),...
    'MarkerFaceColor',colors_detail(col_use,:),'MarkerEdgeColor',colors_main(col_use,:),'LineWidth',1,...
    'CapSize',0);

errorbar(grad_cens,to_mean-to_offset,to_se,'o','MarkerSize',8,'Color',colors_main(1,:),...
    'MarkerFaceColor',colors_detail(1,:),'MarkerEdgeColor',colors_main(1,:),'LineWidth',4,...
    'CapSize',0)

% errorbar(grad_cens,to_mean-to_offset,to_se,'o','MarkerSize',8,'Color',colors_main(1,:),...
%     'MarkerFaceColor',colors_detail(1,:),'MarkerEdgeColor',colors_main(1,:),'LineWidth',4,...
%     'CapSize',0)
% legend([all_pts,wmfit,piplot,ciplot,wpiplot,wciplot],...
%     {'Data','Weighted fit','PI (Unweighted)','CI (unweighted)','PI (weighted)','CI (Weighted)'})


% 
% plot(grads(outlier_mask),to_vals_unc(outlier_mask),'k.','MarkerSize',15)
% errorbar(grad_cens,to_se_cull,...
%     0*to_se_cull,0*to_se_cull,...
%     0.5*diff(grad_edges),0.5*diff(grad_edges),...
%     'r.','MarkerSize',25,'CapSize',0)

xlabel('Gradient of Signal (Hz$^2$/GHz)')
ylabel(sprintf('$f_{TO}$ - %.f (MHz)',to_offset))
ylim([-1e4,1.2e4])
box on
set(gca,'FontSize',20)
set(gca,'LineWidth',2)


plot(check_grad,to_check-to_offset,'bo','MarkerSize',9,'LineWidth',2)

col_use = 1;

% 
% 
% subplot(2,2,3)
hold on
% plot(grads(outlier_mask),to_vals_unc(outlier_mask),'k.','MarkerSize',15)
% errorbar(grad_cens,to_se_cull,...
%     0*to_se_cull,0*to_se_cull,...
%     0.5*diff(grad_edges),0.5*diff(grad_edges),...
%     'r.','MarkerSize',25,'CapSize',0)

set(gca,'FontSize',20)
box on
set(gca,'FontSize',20)
set(gca,'LineWidth',2)
% 
% subplot(2,2,4)
% hold on
% % histogram(mask_out.residuals./to_vals_unc(outlier_mask),linspace(-.15e3,.15e3,60))
% histogram(mask_out.residuals./to_vals_unc(outlier_mask),linspace(-.5,.5,60))
% % errorbar(grad_cens,to_se_cull,...
% %     0*to_se_cull,0*to_se_cull,...
% %     0.5*diff(grad_edges),0.5*diff(grad_edges),...
% %     'r.','MarkerSize',25,'CapSize',0)
% xlabel('Grad')
% ylabel('Variability')
% legend('F unc','binned std')
% set(gca,'FontSize',20)
% box on
% % set(gca,'Yscale','log')
% set(gca,'LineWidth',2)


% subplot(2,2,4)
% hold on
% plot(grads(outlier_mask),grads(outlier_mask).^grad_pow.*to_vals_unc(outlier_mask),'k.','MarkerSize',15)
% % plot(mask_out.residuals,grads(outlier_mask).^grad_pow.*to_vals_unc(outlier_mask),'k.','MarkerSize',15)
% % plot(mask_out.residuals,to_vals_unc(outlier_mask),'r.','MarkerSize',15)
% % xlabel('Residual')
% xlabel('G')
% ylabel(sprintf('G$^%u$ $\\times$ unc.',grad_pow))
% set(gca,'FontSize',20,'Yscale','lin','Xscale','lin')
% % 
% 
% 
% subplot(2,3,6)
% hold on
% plot(grads(outlier_mask),to_vals_unc(outlier_mask),'k.','MarkerSize',15)
% % errorbar(grad_cens,to_se_cull,...
% %     0*to_se_cull,0*to_se_cull,...
% % %     0.5*diff(grad_edges),0.5*diff(grad_edges),...
% %     'r.','MarkerSize',25,'CapSize',0)
% xlabel('Grad')
% ylabel('Variability')
% legend('F unc','binned std')
% set(gca,'FontSize',20)

% Interesting - it seems that the precision (roughly) tracks the variance
% When using weighted fits, does one usually use the standard deviation or
% the variance? The wiki article states: The Guass-Markov theorem shows
% that the model parameter estimate is a best-linear-unbiased-estimator if
% the weight is equal to the reciprocal of the variance; i.e. yes if we use
% the SE from the measurements, then 1/f_unc is the right thing to do.
% And I can see why this might improve the CI - but why does the PI change?
















% Plot binned data
% %
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
% 
% % 
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
% % %
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
% % end