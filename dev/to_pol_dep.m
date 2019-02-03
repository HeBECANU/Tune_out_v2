%quick plot of to dependence on polarisation angle
%hwp_ang = [3,342,353,14,319,330,342,33,37,34,19,320,336,19,313,301,288,20,5,358.5];
%

to_vals_lin_quad = [
    46.5,725735103.8,85,725735152.0,118;
    4.0,725735858.4,54,725735858.2,68;
    27.0,725735290.9,21,725735260.0,24;
    68,725735076.9,35,725735071.0,45;
    137,725735901.3,56,725735962.4,73;
    162, 725735945.6,56,725735944.4,84;
    184.0,725735740.3,53,725735719.0,70;
    103,725735096.7,94,725735106.7,130;
    116,725734919.1,32,725734936.3,57;
    105,725735021.4,77,725734965.0,108;
    75,725734730.4,67,725734879.0,85;
    151,725735811.5,55,725735808.1,73;
    172,725735763.9,72,725735761.9,99;
    75,725735075.2,32,725735059.1,53;
    127,725735937.7,23,725735977.0,29
    101,725735841.7,34,725735857.5,66;
    75,725735479.4,41,725735493.6,71;
    81,725734901.5,52,725734884.5,68;
    51,725735486.3,55,725735493.3,99;
    37,725735574.4,61,725735467.4,82;
    20.5,725735650.9,59,725735645.8,73;
    0,725736176.6,71,725736167.3,67;
    170,725735814.8,60,725735855.9,88;
    150,725736464.8,90,725736482.3,105;
    130.0,725735400.0,77,725735429.8,96;
    107,725735807.8,92,725735642.5,96;
    89,725735587.1,137,725735310.4,118;
    ];%91,725735667.6,15;


%p_min(uw), phi_min(deg) , p_max(uw), phi_max, p_min_perp, phi_min_perp
polz_data = [
    1.87,46.5,257,139,254,135,3.0;
    0.38,4,204,97,169,94,342.0;
    1.0,27,190,116,NaN,118,353.0;
    0.8,68,170,154,170,156,14.0;
    0.58,137.5,167,47,169.7,46,319.0;
    0.15,162,236,72,237,72,330;
    0.5,184,201,273,198,273.5,342.0;
    0.1,103,211,18,203,15,33;
    0.1,116,210,25,206,24,39;
    0.17,105,160,14,NaN,NaN,34;
    0.51,75,150,164,148,165.5,19;
    0.32,151,157,242,241,156,352;
    0.13,172,126,262,119,262,336;
    0.51,75,150,164,148,165.5,19;
    0.92,127,162,217.5,161,218,313;
    0.9,101,143,193,NaN,NaN,301;
    0.14,75,123,164,125,166,288;
    0.38,81,104,170,172,103,20;
    1.03,51,150,137,NaN,NaN,5;
    0.8,37,125.5,123.5,123,124.5,358.5;
    0.72,20.5,123,113,nan,nan,351;
    0.31,0,150,88,150,88,340;
    0.16,170,114,80,116,80,335;
    0.37,150,137,58,nan,nan,324;
    0.76,130,145,39,nan,nan,314;
    0.9,107,153,12.5,nan,nan,303;
    0.3,89,115,0,nan,nan,295];
hwp_ang= polz_data(:,7);
%hwp_ang = hwp_ang-360.*(hwp_ang>180);
%0.32,91,230,7,240,0;
polz_power_frac=2*polz_data(:,1)./(polz_data(:,1)+polz_data(:,3));
%polz_power_frac=2*polz_data(:,1).*polz_data(:,3)./(polz_data(:,1).^2+polz_data(:,3).^2);

bc_angles_wraped=to_vals_lin_quad(:,1);
bc_angles_wraped=bc_angles_wraped+180*(bc_angles_wraped<81)

polz_sign=[0,0,0,0,1,0,0,1,1,1,0,1,0,0,1,1,0,1,0,0,0,0,1,1,1,1,1]';
polz_sign(~polz_sign)=-1;
vec_amp=20e3;

%fit sin waves to the two sets of data
modelfun = @(b,x) b(1).*(cos(x(:,1).*pi/180+b(2).*2*pi).^2)+b(3);
opts = statset('nlinfit');
vec_corr_to=to_vals_lin_quad(:,2)-polz_power_frac.*polz_sign*vec_amp;
ci_size_cut_outliers=1-erf(2.5/sqrt(2));%1-erf(zvalue/sqrt(2)) %confidence interval for cutting outliers
beta0 = [1e3,0.5,nanmean(vec_corr_to)]; %intial guesses
wlin=1./(to_vals_lin_quad(:,3).^2);
fit_mdl_lin = fitnlm(to_vals_lin_quad(:,1),vec_corr_to,modelfun,beta0,...
    'Options',opts,'Weights',wlin,'CoefficientNames' ,{'amp','phase','offset'});

%cut outliers
[~,yci_cull_lim]=predict(fit_mdl_lin,to_vals_lin_quad(:,1),'Prediction','observation','Alpha',ci_size_cut_outliers);
is_outlier_idx=vec_corr_to>yci_cull_lim(:,1) & vec_corr_to<yci_cull_lim(:,2);
vec_corr_to_trim = vec_corr_to(is_outlier_idx);
to_vals_lin_trim = to_vals_lin_quad(is_outlier_idx);

beta0 = fit_mdl_lin.Coefficients{1:3,1}'; %intial guesses
wlin=1./(to_vals_lin_quad(is_outlier_idx,3).^2);
fit_mdl_lin = fitnlm(to_vals_lin_trim,vec_corr_to_trim,modelfun,beta0,...
    'Options',opts,'Weights',wlin,'CoefficientNames' ,{'amp','phase','offset'});
lin_fit_max=[fit_mdl_lin.Coefficients.Estimate(1)+fit_mdl_lin.Coefficients.Estimate(3)...
    sqrt(fit_mdl_lin.Coefficients.SE(1)^2+fit_mdl_lin.Coefficients.SE(3)^2)];

sfigure(866);
clf
ci_size_disp=1-erf(1/sqrt(2));%one sd %confidence interval to display
xsamp = linspace(-50,180,1e4).';
b1 = fit_mdl_lin.Coefficients{1:3,1};
%subplot(2,1,1)
hold on
title('Lin')
xlabel('Input Pol angle (degrees)')
ylabel(sprintf('Tune out value-%.1f±%.1f (MHz)',lin_fit_max(1),lin_fit_max(2)) )
errorbar(to_vals_lin_trim,vec_corr_to_trim-lin_fit_max(1),to_vals_lin_quad(is_outlier_idx,3),'ko')
[y_lin,yci_lin]=predict(fit_mdl_lin,xsamp,'Prediction','observation','Alpha',ci_size_disp); %'Prediction','observation'
%y_lin = b1(1).*(cos(xsamp.*pi/180+b1(2).*2*pi).^2)+b1(3);
plot(xsamp,y_lin-lin_fit_max(1),'b-','LineWidth',1.6)
plot(xsamp,yci_lin(:,1)-lin_fit_max(1),'r-','LineWidth',1.6)
plot(xsamp,yci_lin(:,2)-lin_fit_max(1),'r-','LineWidth',1.6)
scatter(to_vals_lin_quad(~is_outlier_idx),vec_corr_to(~is_outlier_idx)-lin_fit_max(1),60,'rx')

set(gcf,'color','w')
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1600, 900])

% wquad=1./(to_vals_lin(:,3).^2);
% fit_mdl_quad = fitnlm([to_vals_quad(:,1),polz_power_frac],to_vals_quad(:,2),modelfun,beta0,'Options',opts,'Weights',wquad); 
% quad_fit_max=[fit_mdl_quad.Coefficients.Estimate(1)+fit_mdl_quad.Coefficients.Estimate(3)...
%     sqrt(fit_mdl_quad.Coefficients.SE(1)^2+fit_mdl_quad.Coefficients.SE(3)^2)];
% subplot(2,1,2)
% hold on
% b2 = fit_mdl_quad.Coefficients{1:3,1};
% title('Quad')
% xlabel('Input Pol angle (degrees)')
% ylabel(sprintf('Tune out value-%.1f±%.1f (MHz)',quad_fit_max(1),quad_fit_max(2)) )
% errorbar(to_vals_lin(:,1),vec_corr_to-quad_fit_max(1),to_vals_lin(:,3),'ko')
% %[y_quad,yci_quad]=predict(fit_mdl_quad,xsamp,'Prediction','observation','Alpha',ci_size_disp); %'Prediction','observation'
% y_quad = b2(1).*(cos(xsamp.*pi/180+b2(2).*2*pi).^2)+b2(3);
% plot(xsamp,y_lin-quad_fit_max(1),'b-','LineWidth',1.6)
% plot(xsamp,yci_lin(:,1)-quad_fit_max(1),'r-','LineWidth',1.6)
% plot(xsamp,yci_lin(:,2)-quad_fit_max(1),'r-','LineWidth',1.6)

%%
hilight_cut=size(bc_angles_wraped,1)-7;
sin_f = @(b,x) b(1).*(sin(x(:,1).*pi/180+b(2).*2*pi))+b(3);
beta0=[0.015,0,0];
fit_mdl_sin = fitnlm(bc_angles_wraped,-polz_power_frac.*polz_sign,sin_f,beta0,...
    'Options',opts,'CoefficientNames' ,{'amp','phase','offset'});
xsamp = linspace(50,280,1e4).';
[y_lin,yci_lin]=predict(fit_mdl_sin,xsamp,'Prediction','observation','Alpha',ci_size_disp); %'Prediction','observation'
sfigure(553);
clf
scatter(bc_angles_wraped,polz_power_frac,'kx')
xlabel('pol angle')
ylabel('polz power frac')
hold on
scatter(bc_angles_wraped,-polz_power_frac.*polz_sign,'bo')
plot(xsamp,y_lin,'b-','LineWidth',1.6)
plot(xsamp,yci_lin(:,1),'r-','LineWidth',1.6)
plot(xsamp,yci_lin(:,2),'r-','LineWidth',1.6)
set(gcf,'color','w')
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1600, 900])

sfigure(4882);
clf
errorbar(bc_angles_wraped(1:hilight_cut),to_vals_lin_quad(1:hilight_cut,2)-lin_fit_max(1),to_vals_lin_quad(1:hilight_cut,3),'ko')
hold on
errorbar(bc_angles_wraped(hilight_cut+1:end),to_vals_lin_quad(hilight_cut+1:end,2)-lin_fit_max(1),to_vals_lin_quad(hilight_cut+1:end,3),'bx')
xlabel('pol angle')
ylabel('to value')
set(gcf,'color','w')
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1600, 900])

sfigure(6219);
clf
errorbar(hwp_ang(is_outlier_idx),vec_corr_to(is_outlier_idx)-lin_fit_max(1),to_vals_lin_quad(is_outlier_idx,3),'ko')
hold on
errorbar(hwp_ang(~is_outlier_idx),vec_corr_to(~is_outlier_idx)-lin_fit_max(1),to_vals_lin_quad(~is_outlier_idx,3),'rx')
xlabel('hwp angle')
ylabel('to value')
set(gcf,'color','w')
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1600, 900])
