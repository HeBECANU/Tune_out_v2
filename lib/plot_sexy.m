function plot_sexy(disp_config,x_dat,y_dat,weights,fit_mdl)

% Plots "Spectacular Errorbar XY" by binning data in X and producing errorbars.
% Inputs: [X,Y,Y_err] pairs If Y_err empty, defaults so zero
% Configs:
% disp_config.num_bin;
% disp_config.colors_main;
% disp_config.plot_title;
% disp_config.font_name;
% disp_config.font_size_global;
% disp_config.beta0;
% disp_config.opts;
% disp_config.fig_number;

%Bryce comments
% improved integration with code by removing fit functonality

if ~isfield(disp_config,'bin_tol')
    fprintf('no tolerance specified using default\n')
    disp_config.bin_tol=range(x_dat)/1000;
end

colors_main=disp_config.colors_main;
plot_title=disp_config.plot_title;
font_name=disp_config.font_name;
font_size_global=disp_config.font_size_global;

fig_number=disp_config.fig_number;

%% Get data

if isempty(weights)
    weights = 0*y_dat;
end
weights = weights/sum(weights);

if ~isfield(disp_config,'plot_offset') || (~isstruct(disp_config.plot_offset) && disp_config.plot_offset=='avg')
    plot_offset.val=sum(y_dat.*weights)./sum(weights);%predict(fit_mdl,0);
end 
if isstruct(disp_config.plot_offset)
    plot_offset=disp_config.plot_offset;
end


%% Bin up the data
%improved method using uniquetol
bin_center=uniquetol(x_dat,disp_config.bin_tol);
bin_edges=[-inf;(bin_center(2:end)+bin_center(1:end-1))/2;inf];
num_bin=numel(bin_edges)-1;
x_grouped = ones(1,num_bin);
y_grouped=nan(numel(num_bin),3);
x_lims = ones(2,num_bin);
for jj=1:num_bin
    mask=x_dat>bin_edges(jj) & x_dat<bin_edges(jj+1);
    x_grouped(jj) = nanmean(x_dat(mask));
    y_grouped(jj,1)=sum(y_dat(mask).*weights(mask))./sum(weights(mask));
    y_grouped(jj,2)=sqrt(nanvar(y_dat(mask),weights(mask)));
    y_grouped(jj,3)=sewm(y_dat(mask),weights(mask)/sum(weights(mask)));
    try
        x_lims(:,jj) = [min(x_dat(mask)); max(x_dat(mask))];
    catch
        x_lims(:,jj) = [0;0];
    end
end

y_sd = (y_grouped(:,2));
y_se = (y_grouped(:,3));
xneg = x_grouped - x_lims(1,:);
xpos = x_lims(2,:) - x_grouped;

if numel(x_grouped)<2
    error('only a single point was grouped')
end

%% set up plot formatting

colors_main=colors_main./255;
lch=colorspace('RGB->LCH',colors_main(:,:));
lch(:,1)=lch(:,1)+20;
colors_detail=colorspace('LCH->RGB',lch);
%would prefer to use srgb_2_Jab here
color_shaded=colorspace('RGB->LCH',colors_main(3,:));
color_shaded(1)=100;
color_shaded=colorspace('LCH->RGB',color_shaded);



%% Plot the thing!

sfigure(fig_number);
clf
padd_size=0.05*range(x_grouped);
xlim_with_padd=[min(x_grouped)-padd_size,max(x_grouped)+padd_size];
x_samp_pad = linspace(xlim_with_padd(1),xlim_with_padd(2),1e3);
ci_size_disp = 1-erf(1/sqrt(2));
[y_samp_val,y_samp_ci]=predict(fit_mdl,x_samp_pad','Alpha',ci_size_disp); %'Prediction','observation'
%we add another offset so that the plot is about the TO

patch([x_samp_pad, fliplr(x_samp_pad)], ([y_samp_ci(:,1)', fliplr(y_samp_ci(:,2)')]-plot_offset.val), color_shaded,'EdgeColor','none');  %
hold on
%plot(x_samp_pad,(y_samp_ci'-plot_offset.val),'r','color',colors_main(3,:),'LineWidth',1.5);
xl=xlim;
%line(xl,[0,0],'color','k','LineWidth',1)
title(plot_title)
plot(x_samp_pad,(y_samp_val-plot_offset.val),'-','color',colors_main(3,:),'LineWidth',1.5)


errorbar(x_grouped,(y_grouped(:,1)-plot_offset.val),y_sd...
     ,'o','CapSize',0,'Marker','none','Color',colors_detail(1,:),...
     'LineWidth',1.5);
errorbar(x_grouped,(y_grouped(:,1)-plot_offset.val),[],[],...
      xneg,xpos,'o','CapSize',0,'Marker','none','Color',colors_detail(1,:),...
      'LineWidth',1.5);
errorbar(x_grouped,(y_grouped(:,1)-plot_offset.val),y_se,...
    'o','CapSize',0,'MarkerSize',5,'Color',colors_main(1,:),...
    'MarkerFaceColor',colors_detail(1,:),'LineWidth',2.5);
if isfield(plot_offset,'unc')
    ylabel(sprintf('Tune-out value - %.0f\\pm%.0f (MHz)',plot_offset.val,plot_offset.unc))
else
    ylabel(sprintf('Tune-out value - %.3f (MHz)',plot_offset.val))
end

if isfield(disp_config,'x_ticks')
    xticks(disp_config.x_ticks)
end

xlabel(disp_config.x_label)
set(gca,'xlim',xlim_with_padd)
%set(gca,'ylim',first_plot_lims(2,:))
%set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1600, 900])
set(gcf,'color','w')
set(gca,'FontSize',font_size_global,'FontName',font_name)
set(gca, 'Layer','top')
box on

end