function t = corr_plot(x,y,w)
w=w/sum(w);

%Finally plot a nice version of the quad fit
%set up the colors to use
colors_main=[[233,87,0];[33,188,44];[0,165,166]];

colors_main=colors_main./255;
lch=colorspace('RGB->LCH',colors_main(:,:));
lch(:,1)=lch(:,1)+20;
colors_detail=colorspace('LCH->RGB',lch);
%would prefer to use srgb_2_Jab here
color_shaded=colorspace('RGB->LCH',colors_main(3,:));
color_shaded(1)=125;
color_shaded=colorspace('LCH->RGB',color_shaded);

mdl_fun = @(b,x) b(1)+b(2).*x(:,1);
beta0 = [1e14,1e5];
opts = statset('nlinfit');

fit_mdl = fitnlm(x,y,mdl_fun,beta0,'Options',opts,'Weight',w);%'ErrorModel','combined'
ci_size_disp = 1-erf(5/sqrt(2));

x_samp = linspace(min(x), max(x));
[ysamp_culled,yci_culled]=predict(fit_mdl,x_samp','Alpha',ci_size_disp); %'Prediction','observation'
patch([x_samp, fliplr(x_samp)], ([yci_culled(:,1)', fliplr(yci_culled(:,2)')]), color_shaded,'EdgeColor','none');  %
hold on
plot(x_samp,(yci_culled'),'r','color',colors_main(3,:),'LineWidth',1.5);
plot(x_samp,(ysamp_culled),'-','color',colors_main(2,:),'LineWidth',1.5)
scatter(x,y,'kx')
set(gcf,'color','w')
box on
end