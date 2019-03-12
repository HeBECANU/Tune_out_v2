%quick plot of to dependence on polarisation angle
% To predict your new centre freq:
% fprintf('%f',predict(fit_mdl_lin,ANGLE_GOES_HERE,'Prediction','observation','Alpha',ci_size_disp))




to_vals_lin_quad = [
     150,725742096.2,39,725742099.5,51;
     146,725742281.6,25,725742315.2,49;
     142,725742226.3,67,725742125.6,93;
     138,725742346.8,137,725742413.1,153;
     134,725741882.5,108,725742140.2,102;
     130,725741227.9,73,725741142.7,61;
     154,725741961.9,74,725741865.2,91;
     162,725741108.2,94,725740945.7,67;
     177,725738355.7,65,725738400.3,88;
     187,725736236.1,89,725736240.1,121;
     202,725733284.6,76,725733284.0,100;
     220,725730129.5,69,725730014.2,87;
     226,725729582.1,77,725729563.4,106;
     234,725729042.5,19,725728945.3,30;
     246,725728713.1,60,725728802.2,75;
     254,725730171.3,91,725730299.0,154;
     260,725731285.4,53,725731289.0,88;
     286,725736861.9,90,725736679.4,129; %should not be use in final value as has no analog logs
     310,725741328.1,87,725741307.8,104;
     286.5,725737277.3,79,725737412.3,105;
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
     270,725733518.4,83,725733480.3,96
=======
=======
>>>>>>> 6de15477018eba4f08a344e5fafd952c48dbe4ed
=======
>>>>>>> 6de15477018eba4f08a344e5fafd952c48dbe4ed
=======
>>>>>>> 6de15477018eba4f08a344e5fafd952c48dbe4ed
=======
>>>>>>> 6de15477018eba4f08a344e5fafd952c48dbe4ed
     270,725733518.4,83,725733480.3,96;
     246,725729215.0,103,725729418,121;
     283,725735078.9,45,725735160.3,81;
     283,725736295.0,22,725736324.0,35;
     280,725735737.3,36,725735760.8,54
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
>>>>>>> 6de15477018eba4f08a344e5fafd952c48dbe4ed
=======
>>>>>>> 6de15477018eba4f08a344e5fafd952c48dbe4ed
=======
>>>>>>> 6de15477018eba4f08a344e5fafd952c48dbe4ed
=======
>>>>>>> 6de15477018eba4f08a344e5fafd952c48dbe4ed
=======
>>>>>>> 6de15477018eba4f08a344e5fafd952c48dbe4ed
     ];

polz_data = [
   150,90.6,245,63.9,152;
   146,118,260,74,344;
   142,77.5,278,53,9;
   138,134,284,58,3;
   134,156,292,59,18;
   130,119,284,31,12;
   154,113,236,58,328;
   162,132,242,35,329;
   177,150,248,8.7,338;
   187,165,253,0.46,343;
   202,178,262,11.5,355;
   220,110,273,44.1,192;
   226,87,290,73,199;
   234,95,226,82,319;
   246,101,50,42,316;
   254,115,240,22.1,143;
   260,124,51,18.1,329;
   286,1.2,201,0.05,352; %again don't use in final value
   310,140.7,117,40.8,205;
   286.5,165,88.5,2.18,179;
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
   270,167.3,248,5.37,159
=======
=======
>>>>>>> 6de15477018eba4f08a344e5fafd952c48dbe4ed
=======
>>>>>>> 6de15477018eba4f08a344e5fafd952c48dbe4ed
=======
>>>>>>> 6de15477018eba4f08a344e5fafd952c48dbe4ed
=======
>>>>>>> 6de15477018eba4f08a344e5fafd952c48dbe4ed
   270,167.3,248,5.37,159;
   246,127,29,45,136;
   283,195,83,0.6,352;
   283,195,83,0.6,352;
   280,168,257,0.13,171
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
>>>>>>> 6de15477018eba4f08a344e5fafd952c48dbe4ed
=======
>>>>>>> 6de15477018eba4f08a344e5fafd952c48dbe4ed
=======
>>>>>>> 6de15477018eba4f08a344e5fafd952c48dbe4ed
=======
>>>>>>> 6de15477018eba4f08a344e5fafd952c48dbe4ed
=======
>>>>>>> 6de15477018eba4f08a344e5fafd952c48dbe4ed
   ];%

to_lin_val=725736124;


qwp_cen_angle=195;
qwp_ang= polz_data(:,1);
%hwp_ang = hwp_ang-360.*(hwp_ang>180);
%0.32,91,230,7,240,0;
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
A=2*polz_data(:,2).*polz_data(:,4)./(polz_data(:,2).^2+polz_data(:,4).^2);

=======
%A=2*polz_data(:,2).*polz_data(:,4)./(polz_data(:,2).^2+polz_data(:,4).^2);
pol_power_dif=(-polz_data(:,4)+polz_data(:,1))./(polz_data(:,4)+polz_data(:,1));
>>>>>>> 6de15477018eba4f08a344e5fafd952c48dbe4ed
=======
%A=2*polz_data(:,2).*polz_data(:,4)./(polz_data(:,2).^2+polz_data(:,4).^2);
pol_power_dif=(-polz_data(:,4)+polz_data(:,1))./(polz_data(:,4)+polz_data(:,1));
>>>>>>> 6de15477018eba4f08a344e5fafd952c48dbe4ed
=======
%A=2*polz_data(:,2).*polz_data(:,4)./(polz_data(:,2).^2+polz_data(:,4).^2);
pol_power_dif=(-polz_data(:,4)+polz_data(:,1))./(polz_data(:,4)+polz_data(:,1));
>>>>>>> 6de15477018eba4f08a344e5fafd952c48dbe4ed
=======
%A=2*polz_data(:,2).*polz_data(:,4)./(polz_data(:,2).^2+polz_data(:,4).^2);
pol_power_dif=(-polz_data(:,4)+polz_data(:,1))./(polz_data(:,4)+polz_data(:,1));
>>>>>>> 6de15477018eba4f08a344e5fafd952c48dbe4ed
=======
%A=2*polz_data(:,2).*polz_data(:,4)./(polz_data(:,2).^2+polz_data(:,4).^2);
pol_power_dif=(-polz_data(:,4)+polz_data(:,1))./(polz_data(:,4)+polz_data(:,1));
>>>>>>> 6de15477018eba4f08a344e5fafd952c48dbe4ed
A=2*sqrt(polz_data(:,2).*polz_data(:,4))./(polz_data(:,2)+polz_data(:,4));

%%



%fit sin waves to the two sets of data
modelfun = @(b,x) b(1).*sin(x(:,1).*2*pi./180 -b(3).*(2*pi/360) )+b(2);
opts = statset('MaxIter',1e4);
beta0 = [6e3,to_lin_val,qwp_cen_angle]; %intial guesses
wlin=1./(to_vals_lin_quad(:,3).^2);
fit_mdl_lin = fitnlm(qwp_ang,to_vals_lin_quad(:,2),modelfun,beta0,...
    'Options',opts,'Weights',wlin,'CoefficientNames' ,{'amp','freq_offset','fast_angle'});


colors_main=[[233,87,0];[33,188,44];[0,165,166]];
colors_main=colors_main./255;
lch=colorspace('RGB->LCH',colors_main(:,:));
lch(:,1)=lch(:,1)+20;
colors_detail=colorspace('LCH->RGB',lch);
%would prefer to use srgb_2_Jab here
color_shaded=colorspace('RGB->LCH',colors_main(3,:));
color_shaded(1)=100;
color_shaded=colorspace('LCH->RGB',color_shaded);

<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD



=======
>>>>>>> 6de15477018eba4f08a344e5fafd952c48dbe4ed
=======
>>>>>>> 6de15477018eba4f08a344e5fafd952c48dbe4ed
=======
>>>>>>> 6de15477018eba4f08a344e5fafd952c48dbe4ed
=======
>>>>>>> 6de15477018eba4f08a344e5fafd952c48dbe4ed
sfigure(3)
clf
set(gcf,'color','w')
xlabel('QWP angle (º)')
ylabel(sprintf('Tune out value Relative to Linear Pol (MHz)'))
colors_main=[[233,87,0];[33,188,44];[0,165,166]];
colors_main=colors_main./255;
lch=colorspace('RGB->LCH',colors_main(:,:));
lch(:,1)=lch(:,1)+20;
colors_detail=colorspace('LCH->RGB',lch);
%would prefer to use srgb_2_Jab here
color_shaded=colorspace('RGB->LCH',colors_main(3,:));
color_shaded(1)=100;
color_shaded=colorspace('LCH->RGB',color_shaded);
ci_size_disp=1-erf(1/sqrt(2));%one sd %confidence interval to display
plot_padd=50;
xsamp = linspace(min(qwp_ang)-plot_padd,max(qwp_ang)+plot_padd,1e4).';
fit_values = fit_mdl_lin.Coefficients{:,1};
%subplot(2,1,1)
hold on
[y_lin,yci_lin]=predict(fit_mdl_lin,xsamp,'Prediction','observation','Alpha',ci_size_disp); %'Prediction','observation'
%y_lin = b1(1).*(cos(xsamp.*pi/180+b1(2).*2*pi).^2)+b1(3);
plot(xsamp,y_lin-to_lin_val,'b-','LineWidth',1.6)
plot(xsamp,yci_lin(:,1)-to_lin_val,'-','LineWidth',1.6,'Color',[1,1,1]*0.6)
plot(xsamp,yci_lin(:,2)-to_lin_val,'-','LineWidth',1.6,'Color',[1,1,1]*0.6)
errorbar(qwp_ang,to_vals_lin_quad(:,2)-to_lin_val,to_vals_lin_quad(:,3),'o','CapSize',0,'MarkerSize',5,'Color',colors_main(1,:),...
    'MarkerFaceColor',colors_detail(1,:),'LineWidth',2.5)
hold off
box on
xlim([120,320])
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1600, 900])


[~,A_min_idx] = min(A);
counts = 1:length(A);
sign_flip = ones(size(A')).*(1-2*(counts>A_min_idx));
sign_flip = ones(size(A')).*(1-2*(45>abs(polz_data(:,1)'-235)));
A_plot = A.*sign_flip';

sfigure(7920);
clf;
subplot(2,2,1)
scatter(polz_data(:,1),A_plot)
title('A vz angle')
subplot(2,2,2)
scatter(polz_data(:,1),to_vals_lin_quad(:,2))
title('TO vs angle')
subplot(2,2,3)
scatter(A_plot,to_vals_lin_quad(:,2))
title('TO vs A')
suptitle('TO dependence on ?pplz')
subplot(2,2,4)
scatter(polz_data(:,1),unwrap(mod(polz_data(:,5)+45,180).*pi/180))
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
title('qwp vs min pol ang')
=======
=======
>>>>>>> 6de15477018eba4f08a344e5fafd952c48dbe4ed
=======
>>>>>>> 6de15477018eba4f08a344e5fafd952c48dbe4ed
=======
>>>>>>> 6de15477018eba4f08a344e5fafd952c48dbe4ed
=======
>>>>>>> 6de15477018eba4f08a344e5fafd952c48dbe4ed
title('qwp vs min pol ang')
sfigure(2782);
subplot(2,1,1)
scatter(polz_data(:,1),pol_power_dif)
xlabel('QWP ang')
ylabel('pol power dif')
subplot(2,1,2)
scatter(polz_data(:,1),mod(polz_data(:,5),180))
xlabel('QWP ang')
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
ylabel('min power ang')
>>>>>>> 6de15477018eba4f08a344e5fafd952c48dbe4ed
=======
ylabel('min power ang')
>>>>>>> 6de15477018eba4f08a344e5fafd952c48dbe4ed
=======
ylabel('min power ang')
>>>>>>> 6de15477018eba4f08a344e5fafd952c48dbe4ed
=======
ylabel('min power ang')
>>>>>>> 6de15477018eba4f08a344e5fafd952c48dbe4ed
=======
ylabel('min power ang')
>>>>>>> 6de15477018eba4f08a344e5fafd952c48dbe4ed
