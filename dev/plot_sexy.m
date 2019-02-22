function plot_sexy(disp_config,X,Y,Y_err)

% Plots "Spectacular Errorbar XY" by binning data in X and producing errorbars.
% Inputs: [X,Y,Y_err] pairs If Y_err empty, defaults so zero
% Configs:
% disp_config.num_bin;
% disp_config.colors_main;
% disp_config.plot_title;
% disp_config.font_name;
% disp_config.font_size_global;
% disp_config.mdl_fun;
% disp_config.beta0;
% disp_config.opts;
% disp_config.fig_number;

num_bin=disp_config.num_bin;

colors_main=disp_config.colors_main;
plot_title=disp_config.plot_title;
font_name=disp_config.font_name;
font_size_global=disp_config.font_size_global;

mdl_fun=disp_config.mdl_fun;
beta0=disp_config.beta0;
opts=disp_config.opts;
fig_number=disp_config.fig_number;

    %% Get data

    if ~isempty(Y_err)
        Y_err=Y_err;
    else
        Y_err = 0*Y;
    end

%     %Plot the raw data
%     sfigure(6695)
%     errorbar(X,Y,Y_err,'k.','CapSize',0)

    % Set fit weights
    weights = 1./Y.^2;
    weights = weights/sum(weights);


    %% Bin up the data

    x_range = max(X)-min(X);
    bin_size = x_range/num_bin;
    c_data = viridis(num_bin);
    x_grouped = ones(1,num_bin);
    y_grouped=nan(numel(num_bin),2);
    x_lims = ones(2,num_bin);
    for jj=0:(num_bin-1)
        bin_centre = bin_size*0.5+jj*bin_size+min(X);
        x_mask = (abs(X-bin_centre)<(bin_size*0.5));
        x_grouped(jj+1) = nanmean(X(x_mask));
        y_grouped(jj+1,1)=sum(Y(x_mask).*weights(x_mask))./sum(weights(x_mask));
        y_grouped(jj+1,2)=sqrt(nanvar(Y(x_mask),weights(x_mask)));
        try
            x_lims(:,jj+1) = [min(X(x_mask)); max(X(x_mask))];
        catch
            x_lims(:,jj+1) = [0;0];
        end
    end

    yneg = (y_grouped(:,2));
    ypos = (y_grouped(:,2));
    xneg = x_grouped - x_lims(1,:);
    xpos = x_lims(2,:) - x_grouped;


    %% set up plot formatting

    colors_main=colors_main./255;
    lch=colorspace('RGB->LCH',colors_main(:,:));
    lch(:,1)=lch(:,1)+20;
    colors_detail=colorspace('LCH->RGB',lch);
    %would prefer to use srgb_2_Jab here
    color_shaded=colorspace('RGB->LCH',colors_main(3,:));
    color_shaded(1)=125;
    color_shaded=colorspace('LCH->RGB',color_shaded);

    fit_mdl = fitnlm(X',Y',mdl_fun,beta0,'Options',opts,'Weight',weights);%'ErrorModel','combined'
    ci_size_disp = 1-erf(1/sqrt(2));


    %% Plot the thing!

    sfigure(fig_number);
    clf
    x_grouped_pad = linspace(-1, x_grouped(end)*2.01,6);
    [ysamp_culled,yci_culled]=predict(fit_mdl,x_grouped_pad','Alpha',ci_size_disp); %'Prediction','observation'
    %we add another offset so that the plot is about the TO
    plot_offset=sum(Y.*weights)./sum(weights);%predict(fit_mdl,0);
    patch([x_grouped_pad, fliplr(x_grouped_pad)], ([yci_culled(:,1)', fliplr(yci_culled(:,2)')]-plot_offset), color_shaded,'EdgeColor','none');  %
    hold on
    plot(x_grouped_pad,(yci_culled'-plot_offset),'r','color',colors_main(3,:),'LineWidth',1.5);
    xl=xlim;
    %line(xl,[0,0],'color','k','LineWidth',1)
    title(plot_title)
    plot(x_grouped_pad,(ysamp_culled-plot_offset),'-','color',colors_main(2,:),'LineWidth',1.5)

    errorbar(x_grouped,(y_grouped(:,1)-plot_offset),yneg,ypos,xneg,xpos,'o','CapSize',0,'MarkerSize',5,'Color',colors_main(1,:),...
        'MarkerFaceColor',colors_detail(1,:),'LineWidth',1.5);
    xlabel('Amplitude of Oscillations (mm)')
    ylabel(sprintf('Tune-out value - %.3f (MHz)',plot_offset))
    set(gca,'xlim',[1,ceil(max(x_grouped)*1.05)])
    %set(gca,'ylim',first_plot_lims(2,:))
    set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1600, 900])
    set(gcf,'color','w')
    set(gca,'FontSize',font_size_global,'FontName',font_name)
    box on

end