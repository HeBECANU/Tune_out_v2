%quick plot of to dependence on polarisation angle
%hwp_ang = [3,342,353,14,319,330,342,33,37,34,19,320,336,19,313,301,288,20,5,358.5];
%
%all data
% to_vals_lin_quad = [
%     46.5,725735103.8,85,725735152.0,118;
%     4.0,725735858.4,54,725735858.2,68;
%     27.0,725735290.9,21,725735260.0,24;
%     68,725735076.9,35,725735071.0,45;
%     137,725735901.3,56,725735962.4,73;
%     162, 725735945.6,56,725735944.4,84;
%     184.0,725735740.3,53,725735719.0,70;
%     103,725735096.7,94,725735106.7,130;
%     116,725734919.1,32,725734936.3,57;
%     105,725735021.4,77,725734965.0,108;
%     75,725734730.4,67,725734879.0,85;
%     151,725735811.5,55,725735808.1,73;
%     172,725735763.9,72,725735761.9,99;
%     75,725735075.2,32,725735059.1,53;
%     127,725735937.7,23,725735977.0,29
%     101,725735841.7,34,725735857.5,66;
%     75,725735479.4,41,725735493.6,71;
%     81,725734901.5,52,725734884.5,68;
%     51,725735486.3,55,725735493.3,99;
%     37,725735574.4,61,725735467.4,82;
%     20.5,725735650.9,59,725735645.8,73;
%     0,725736176.6,71,725736167.3,67;
%     170,725735814.8,60,725735855.9,88;
%     150,725736464.8,90,725736482.3,105;
%     130.0,725735400.0,77,725735429.8,96;
%     107,725735807.8,92,725735642.5,96;
%     89,725735587.1,137,725735310.4,118;
%     111.0,725735526.5,80,725735593.9,109;
%     130.5,725735904.6,66,725735829.1,62;
%     151.0,725736107.4,36,725736071.3,49;
%     170,725736032.5,56,725736054.5,67;
%     190,725735463.0,71,725735467.0,108;
%     29,725735392.5,24,725735337.3,34;
%     50.5,725735068.4,81,725734979.2,82;
%     70,725734755.9,91,725734885.7,131;
%     92,725734989.4,87,725734999.8,116;
%     61,725735043.9,76,725735060.2,92;
%     46,725735483.9,71,725735473.4,117;
%     24,725735498.6,70,725735408.5,103;
%     180,725736189.9,76,725736160.3,138;
%     160,725736304.2,65,725736338.6,96;
%     140.5,725735936.7,28,725736022.8,33;
%     121,725735912.7,20,725735902.6,19;
%     100,725735794.7,53,725735667.8,86;
%     80,725735283.9,45,725735284.2,68];%91,725735667.6,15;
% 
% %p_min(uw), phi_min(deg) , p_max(uw), phi_max, p_min_perp, phi_min_perp
% polz_data = [
%     1.87,46.5,257,139,254,135,3.0;
%     0.38,4,204,97,169,94,342.0;
%     1.0,27,190,116,NaN,118,353.0;
%     0.8,68,170,154,170,156,14.0;
%     0.58,137.5,167,47,169.7,46,319.0;
%     0.15,162,236,72,237,72,330;
%     0.5,184,201,273,198,273.5,342.0;
%     0.1,103,211,18,203,15,33;
%     0.1,116,210,25,206,24,39;
%     0.17,105,160,14,NaN,NaN,34;
%     0.51,75,150,164,148,165.5,19;
%     0.32,151,157,242,241,156,325;
%     0.13,172,126,262,119,262,336;
%     0.51,75,150,164,148,165.5,19;
%     0.92,127,162,217.5,161,218,313;
%     0.9,101,143,193,NaN,NaN,301;
%     0.14,75,123,164,125,166,288;
%     0.38,81,104,170,172,103,20;
%     1.03,51,150,137,NaN,NaN,5;
%     0.8,37,125.5,123.5,123,124.5,358.5;
%     0.72,20.5,123,113,nan,nan,351;
%     0.31,0,150,88,150,88,340;
%     0.16,170,114,80,116,80,335;
%     0.37,150,137,58,nan,nan,324;
%     0.76,130,145,39,nan,nan,314;
%     0.9,107,153,12.5,nan,nan,303;
%     0.3,89,115,0,nan,nan,295;
%     0.34,111.0,79.0,20.0,NaN,NaN,304;
%     0.37,130.5,120.0,40.0,nan,nan,314;
%     0.14,151.0,98.0,63.0,nan,nan,325.0;
%     0.13,170,113,84,nan,nan,335;
%     0.68,190,133,104,nan,nan,344;
%     0.94,29,89,120,85,119,354;
%     0.96,50.5,75,146,nan,nan,6;
%     0.59,70,120,165,nan,nan,15;
%     0.06,92,77,181,nan,nan,26;
%     0.52,61,84,154,nan,nan,10;
%     0.86,46,74,137,nan,nan,3;
%     0.54,24.5,72,112,nan,nan,352;
%     0.19,180,93,89,nan,nan,339;
%     0.12,160,120,70.5,nan,nan,319;
%     0.33,140.5,112,50,118,50,320;
%     0.36,121,86,32,nan,nan,310;
%     0.17,100,72,191,nan,nan,300;
%     0.07,80,107.5,165,nan,nan,290];

%old data( pre nuller reconfig)
% to_vals_lin_quad = [
%     46.5,725735103.8,85,725735152.0,118;
%     4.0,725735858.4,54,725735858.2,68;
%     27.0,725735290.9,21,725735260.0,24;
%     68,725735076.9,35,725735071.0,45;
%     137,725735901.3,56,725735962.4,73;
%     162, 725735945.6,56,725735944.4,84;
%     184.0,725735740.3,53,725735719.0,70;
%     103,725735096.7,94,725735106.7,130;
%     116,725734919.1,32,725734936.3,57;
%     105,725735021.4,77,725734965.0,108;
%     75,725734730.4,67,725734879.0,85;
%     151,725735811.5,55,725735808.1,73;
%     172,725735763.9,72,725735761.9,99;
%     75,725735075.2,32,725735059.1,53;
%     127,725735937.7,23,725735977.0,29
%     101,725735841.7,34,725735857.5,66;
%     75,725735479.4,41,725735493.6,71;
%     81,725734901.5,52,725734884.5,68;
%     51,725735486.3,55,725735493.3,99;
%     37,725735574.4,61,725735467.4,82;
%     20.5,725735650.9,59,725735645.8,73;
%     0,725736176.6,71,725736167.3,67;
%     170,725735814.8,60,725735855.9,88;
%     150,725736464.8,90,725736482.3,105;
%     130.0,725735400.0,77,725735429.8,96;
%     107,725735807.8,92,725735642.5,96;
%     89,725735587.1,137,725735310.4,118];%91,725735667.6,15;
% 
% %p_min(uw), phi_min(deg) , p_max(uw), phi_max, p_min_perp, phi_min_perp
% polz_data = [
%     1.87,46.5,257,139,254,135,3.0;
%     0.38,4,204,97,169,94,342.0;
%     1.0,27,190,116,NaN,118,353.0;
%     0.8,68,170,154,170,156,14.0;
%     0.58,137.5,167,47,169.7,46,319.0;
%     0.15,162,236,72,237,72,330;
%     0.5,184,201,273,198,273.5,342.0;
%     0.1,103,211,18,203,15,33;
%     0.1,116,210,25,206,24,39;
%     0.17,105,160,14,NaN,NaN,34;
%     0.51,75,150,164,148,165.5,19;
%     0.32,151,157,242,241,156,325;
%     0.13,172,126,262,119,262,336;
%     0.51,75,150,164,148,165.5,19;
%     0.92,127,162,217.5,161,218,313;
%     0.9,101,143,193,NaN,NaN,301;
%     0.14,75,123,164,125,166,288;
%     0.38,81,104,170,172,103,20;
%     1.03,51,150,137,NaN,NaN,5;
%     0.8,37,125.5,123.5,123,124.5,358.5;
%     0.72,20.5,123,113,nan,nan,351;
%     0.31,0,150,88,150,88,340;
%     0.16,170,114,80,116,80,335;
%     0.37,150,137,58,nan,nan,324;
%     0.76,130,145,39,nan,nan,314;
%     0.9,107,153,12.5,nan,nan,303;
%     0.3,89,115,0,nan,nan,295
%     ];

%data since nuller reconfig
to_vals_lin_quad = [
    111.0,725735526.5,80,725735593.9,109;
    130.5,725735904.6,66,725735829.1,62;
    151.0,725736107.4,36,725736071.3,49;
    170,725736032.5,56,725736054.5,67;
    190,725735463.0,71,725735467.0,108;
    29,725735392.5,24,725735337.3,34;
    50.5,725735068.4,81,725734979.2,82;
    70,725734755.9,91,725734885.7,131;
    92,725734989.4,87,725734999.8,116;
    61,725735043.9,76,725735060.2,92;
    46,725735483.9,71,725735473.4,117;
    24,725735498.6,70,725735408.5,103;
    180,725736189.9,76,725736160.3,138;
    160,725736304.2,65,725736338.6,96;
    140.5,725735936.7,28,725736022.8,33;
    121,725735912.7,20,725735902.6,19;
    100,725735794.7,53,725735667.8,86;
    80,725735283.9,45,725735284.2,68;
    99.5,725735478.2,110,725735449.7,183;
    120,725735863.3,96,725735952.4,143;
    140,725736123.5,79,725736002.3,111;
    145,725735989.4,112,725736008.6,195;
    155,725736124.9,83,725736053.9,117;
    160,725736319.1,111,725736317.5,123;
    165,725736168.1,92,725736256.5,83;
    171,725736412.3,110,725736476.3,127;
    176.5,725736333.6,136,725736250.3,186;
    176.5,725736175.1,126,725736102.8,156;
    181,725736038.3,92,725735944.2,134;
    187,725736002.5,106,725735973.1,119;
    194,725735718.1,71,725735824.9,100;
    199,725735955.7,91,725736048.5,76;
    208,725735856.0,121,725735892.4,179;
    217,725735356.3,86,725735488.2,93;
    231,725735378.1,95,725735416.3,130;
    239.5,725735501.2,154,725735800.5,209;
    249,725734741.3,171,725734822.9,214];


polz_data = [
    0.34,111.0,79.0,20.0,NaN,NaN,304;
    0.37,130.5,120.0,40.0,nan,nan,314;
    0.14,151.0,98.0,63.0,nan,nan,325.0;
    0.13,170,113,84,nan,nan,335;
    0.68,190,133,104,nan,nan,344;
    0.94,29,89,120,85,119,354;
    0.96,50.5,75,146,nan,nan,6;
    0.59,70,120,165,nan,nan,15;
    0.06,92,77,181,nan,nan,26;
    0.52,61,84,154,nan,nan,10;
    0.86,46,74,137,nan,nan,3;
    0.54,24.5,72,112,nan,nan,352;
    0.19,180,93,89,nan,nan,339;
    0.12,160,120,70.5,nan,nan,319;
    0.33,140.5,112,50,118,50,320;
    0.36,121,86,32,nan,nan,310;
    0.17,100,72,191,nan,nan,300;
    0.07,80,107.5,165,nan,nan,290;
    0.25,99.5,72,188,72,187.5,299;
    0.37,120,78,221,nan,nan,309;
    0.24,140,62,230,nan,nan,320;
    0.23,147,87.6,238,nan,nan,324;
    0.11,155,60.4,248,nan,nan,327;
    0.14,160,76,249,nan,nan,329;
    0.20,165,80,255,77,254,332;
    0.16,171,79,264,nan,nan,335;
    0.23,176.5,65,266,nan,nan,338.5;
    0.23,176.5,65,266,nan,nan,338.5;
    0.33,181,71,271,nan,nan,340;
    0.52,187,92,278,nan,nan,342;
    0.66,194,78,282,nan,nan,347;
    0.72,199,82,290,80,291,350;
    0.65,208,57,298,nan,nan,353;
    0.90,217,78,308,nan,nan,358;
    0.70,231,70.1,135,nan,nan,4;
    0.64,239.5,53,150,60,151,9.5;
    0.58,249,75.1,160,nan,nan,14.5];


hwp_ang= polz_data(:,7);
%hwp_ang = hwp_ang-360.*(hwp_ang>180);
%0.32,91,230,7,240,0;
polz_power_frac=2*polz_data(:,1)./(polz_data(:,1)+polz_data(:,3));
%polz_power_frac=2*polz_data(:,1).*polz_data(:,3)./(polz_data(:,1).^2+polz_data(:,3).^2);

bc_angles_wraped=polz_data(:,2);
%move the data into the center
shift=-(171-90);
shift=-(171-90);
bc_angles_wraped=mod(bc_angles_wraped+shift,180)-shift;

%polz_sign=[0,0,0,0,1,0,0,1,1,1,0,1,0,0,1,1,0,1,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,1,0,0,0,0,1,1,1,1,0]';
%polz_sign=[0,0,0,0,1,0,0,1,1,1,0,1,0,0,1,1,0,1,0,0,0,0,1,1,1,1,1]';
polz_sign=[1,1,1,0,0,0,0,0,1,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]';
polz_sign(~polz_sign)=-1;
vec_amp=0e3%19e3;

%% FIt and plot (with freq as free)
%fit sin waves to the two sets of data
modelfun = @(b,x) b(1).*(cos(x(:,1).*pi./b(4)+b(2).*2*pi).^2)+b(3);
opts = statset('nlinfit');
vec_corr_to=to_vals_lin_quad(:,2)-polz_power_frac.*polz_sign*vec_amp;
ci_size_cut_outliers=1-erf(10/sqrt(2));%1-erf(zvalue/sqrt(2)) %confidence interval for cutting outliers
beta0 = [1e3,0.5,nanmean(vec_corr_to),180]; %intial guesses
wlin=1./(to_vals_lin_quad(:,3).^2);
fit_mdl_lin = fitnlm(polz_data(:,2),vec_corr_to,modelfun,beta0,...
    'Options',opts,'Weights',wlin,'CoefficientNames' ,{'amp','phase','offset','period'});

%cut outliers
[~,yci_cull_lim]=predict(fit_mdl_lin,to_vals_lin_quad(:,1),'Prediction','observation','Alpha',ci_size_cut_outliers);
is_outlier_idx=vec_corr_to>yci_cull_lim(:,1) & vec_corr_to<yci_cull_lim(:,2);
vec_corr_to_trim = vec_corr_to(is_outlier_idx);
to_vals_lin_trim = to_vals_lin_quad(is_outlier_idx);

beta0 = fit_mdl_lin.Coefficients{1:4,1}'; %intial guesses
wlin=1./(to_vals_lin_quad(is_outlier_idx,3).^2);
fit_mdl_lin = fitnlm(to_vals_lin_trim,vec_corr_to_trim,modelfun,beta0,...
    'Options',opts,'Weights',wlin,'CoefficientNames' ,{'amp','phase','offset','period'});
lin_fit_max=[fit_mdl_lin.Coefficients.Estimate(1)+fit_mdl_lin.Coefficients.Estimate(3)...
    sqrt(fit_mdl_lin.Coefficients.SE(1)^2+fit_mdl_lin.Coefficients.SE(3)^2)];

sfigure(866);
clf
ci_size_disp=1-erf(1/sqrt(2));%one sd %confidence interval to display
xsamp = linspace(-50,220,1e4).';
fit_values = fit_mdl_lin.Coefficients{1:4,1};
%subplot(2,1,1)
hold on
title(['amp= ',num2str(fit_values(1)),'MHz , phase= ',num2str(fit_values(2)),', offset= ',num2str(fit_values(3)),'MHz , period= ',num2str(fit_values(4)),'^\circ'])
xlabel('Input Pol angle (degrees)')
ylabel(sprintf('Tune out value-%.1f±%.1f (MHz)',lin_fit_max(1),lin_fit_max(2)) )
errorbar(to_vals_lin_trim,vec_corr_to_trim-lin_fit_max(1),to_vals_lin_quad(is_outlier_idx,3),'ko')
[y_lin,yci_lin]=predict(fit_mdl_lin,xsamp,'Prediction','observation','Alpha',ci_size_disp); %'Prediction','observation'
%y_lin = b1(1).*(cos(xsamp.*pi/180+b1(2).*2*pi).^2)+b1(3);
plot(xsamp,y_lin-lin_fit_max(1),'b-','LineWidth',1.6)
plot(xsamp,yci_lin(:,1)-lin_fit_max(1),'r-','LineWidth',1.6)
plot(xsamp,yci_lin(:,2)-lin_fit_max(1),'r-','LineWidth',1.6)
errorbar(to_vals_lin_quad(~is_outlier_idx),vec_corr_to(~is_outlier_idx)-lin_fit_max(1),to_vals_lin_quad(~is_outlier_idx,3),'rx')

set(gcf,'color','w')

%%
%% FIt and plot (with freq as free)
%fit sin waves to the two sets of data
modelfun = @(b,x) b(1).*(cos(x(:,1).*pi./b(4)+b(2).*2*pi).^2)+b(3);
opts = statset('nlinfit');
vec_corr_to=to_vals_lin_quad(:,2)-polz_power_frac.*polz_sign*vec_amp;
ci_size_cut_outliers=1-erf(10/sqrt(2));%1-erf(zvalue/sqrt(2)) %confidence interval for cutting outliers
beta0 = [1e3,0.5,nanmean(vec_corr_to),180]; %intial guesses
wlin=1./(to_vals_lin_quad(:,3).^2);
fit_mdl_lin = fitnlm(polz_data(:,2),vec_corr_to,modelfun,beta0,...
    'Options',opts,'Weights',wlin,'CoefficientNames' ,{'amp','phase','offset','period'});

%cut outliers
[~,yci_cull_lim]=predict(fit_mdl_lin,to_vals_lin_quad(:,1),'Prediction','observation','Alpha',ci_size_cut_outliers);
is_outlier_idx=vec_corr_to>yci_cull_lim(:,1) & vec_corr_to<yci_cull_lim(:,2);
vec_corr_to_trim = vec_corr_to(is_outlier_idx);
to_vals_lin_trim = to_vals_lin_quad(is_outlier_idx);

beta0 = fit_mdl_lin.Coefficients{1:4,1}'; %intial guesses
wlin=1./(to_vals_lin_quad(is_outlier_idx,3).^2);
%fit_mdl_lin = fitnlm(to_vals_lin_trim,vec_corr_to_trim,modelfun,beta0,...
%    'Options',opts,'Weights',wlin,'CoefficientNames' ,{'amp','phase','offset','period'});
fit_mdl_lin = fitnlm(to_vals_lin_trim,vec_corr_to_trim,modelfun,beta0,...
    'Options',opts,'CoefficientNames' ,{'amp','phase','offset','period'});
lin_fit_max=[fit_mdl_lin.Coefficients.Estimate(1)+fit_mdl_lin.Coefficients.Estimate(3)...
    sqrt(fit_mdl_lin.Coefficients.SE(1)^2+fit_mdl_lin.Coefficients.SE(3)^2)];

sfigure(866);
clf
ci_size_disp=1-erf(1/sqrt(2));%one sd %confidence interval to display
xsamp = linspace(-50,220,1e4).';
fit_values = fit_mdl_lin.Coefficients{1:4,1};
%subplot(2,1,1)
hold on
title(['amp= ',num2str(fit_values(1)),'MHz , phase= ',num2str(fit_values(2)),', offset= ',num2str(fit_values(3)),'MHz , period= ',num2str(fit_values(4)),'^\circ'])
xlabel('Input Pol angle (degrees)')
ylabel(sprintf('Tune out value-%.1f±%.1f (MHz)',lin_fit_max(1),lin_fit_max(2)) )
errorbar(to_vals_lin_trim,vec_corr_to_trim-lin_fit_max(1),to_vals_lin_quad(is_outlier_idx,3),'ko')
[y_lin,yci_lin]=predict(fit_mdl_lin,xsamp,'Prediction','observation','Alpha',ci_size_disp); %'Prediction','observation'
%y_lin = b1(1).*(cos(xsamp.*pi/180+b1(2).*2*pi).^2)+b1(3);
plot(xsamp,y_lin-lin_fit_max(1),'b-','LineWidth',1.6)
plot(xsamp,yci_lin(:,1)-lin_fit_max(1),'r-','LineWidth',1.6)
plot(xsamp,yci_lin(:,2)-lin_fit_max(1),'r-','LineWidth',1.6)
errorbar(to_vals_lin_quad(~is_outlier_idx),vec_corr_to(~is_outlier_idx)-lin_fit_max(1),to_vals_lin_quad(~is_outlier_idx,3),'rx')

set(gcf,'color','w')

%% FIt and plot (with freq fixed)
fprintf('Fixed period fit\n')

%fit sin waves to the two sets of data
modelfun = @(b,x) b(1).*(cos(x(:,1).*pi./180+b(2).*2*pi).^2)+b(3);
opts = statset('MaxIter',1e4);
ci_size_cut_outliers=1-erf(10/sqrt(2));%1-erf(zvalue/sqrt(2)) %confidence interval for cutting outliers
beta0 = [1e3,0.5,nanmean(vec_corr_to)]; %intial guesses
wlin=1./(to_vals_lin_quad(:,3).^2);
fit_mdl_lin = fitnlm(bc_angles_wraped,vec_corr_to,modelfun,beta0,...
    'Options',opts,'Weights',wlin,'CoefficientNames' ,{'amp','phase','offset'});

%cut outliers
[~,yci_cull_lim]=predict(fit_mdl_lin,bc_angles_wraped,'Prediction','observation','Alpha',ci_size_cut_outliers);
is_outlier_idx=vec_corr_to>yci_cull_lim(:,1) & vec_corr_to<yci_cull_lim(:,2);
vec_corr_to_trim = vec_corr_to(is_outlier_idx);
to_vals_lin_trim = to_vals_lin_quad(is_outlier_idx);

beta0 = fit_mdl_lin.Coefficients{:,1}'; %intial guesses
%wlin=1./(to_vals_lin_quad(is_outlier_idx,3).^2);
% fit_mdl_lin = fitnlm(bc_angles_wraped(is_outlier_idx),vec_corr_to(is_outlier_idx),modelfun,beta0,...
%     'Options',opts,'Weights',wlin,'CoefficientNames' ,{'amp','phase','offset'});
fit_mdl_lin = fitnlm(bc_angles_wraped(is_outlier_idx),vec_corr_to(is_outlier_idx),modelfun,beta0,...
    'Options',opts,'CoefficientNames' ,{'amp','phase','offset'});
lin_fit_max=[fit_mdl_lin.Coefficients.Estimate(1)+fit_mdl_lin.Coefficients.Estimate(3)...
    sqrt(fit_mdl_lin.Coefficients.SE(1)^2+fit_mdl_lin.Coefficients.SE(3)^2)];

sfigure(867);
clf
ci_size_disp=1-erf(1/sqrt(2));%one sd %confidence interval to display
plot_padd=20;
xsamp = linspace(min(bc_angles_wraped)-plot_padd,max(bc_angles_wraped)+plot_padd,1e4).';
fit_values = fit_mdl_lin.Coefficients{:,1};
%subplot(2,1,1)
hold on
title(['amp= ',num2str(fit_values(1)),'MHz , phase= ',num2str(fit_values(2)),', offset= ',num2str(fit_values(3)),'MHz , period=180(fixed)^\circ'])
xlabel('Input Pol angle (degrees)')
ylabel(sprintf('Tune out value-%.1f±%.1f (MHz)',lin_fit_max(1),lin_fit_max(2)) )
fprintf('Tune out value-%.1f±%.1f (MHz) (BLUE)\n',lin_fit_max(1),lin_fit_max(2))
fprintf('Tune out value-%.1f±%.1f (MHz) (RED)\n',lin_fit_max(1)/2,lin_fit_max(2)/2)
errorbar(bc_angles_wraped,vec_corr_to_trim-lin_fit_max(1),to_vals_lin_quad(is_outlier_idx,3),'ko')
[y_lin,yci_lin]=predict(fit_mdl_lin,xsamp,'Prediction','observation','Alpha',ci_size_disp); %'Prediction','observation'
%y_lin = b1(1).*(cos(xsamp.*pi/180+b1(2).*2*pi).^2)+b1(3);
plot(xsamp,y_lin-lin_fit_max(1),'b-','LineWidth',1.6)
plot(xsamp,yci_lin(:,1)-lin_fit_max(1),'r-','LineWidth',1.6)
plot(xsamp,yci_lin(:,2)-lin_fit_max(1),'r-','LineWidth',1.6)
errorbar(to_vals_lin_quad(~is_outlier_idx),vec_corr_to(~is_outlier_idx)-lin_fit_max(1),to_vals_lin_quad(~is_outlier_idx,3),'rx')

set(gcf,'color','w')

%% Plot wrapped data
hilight_cut=size(bc_angles_wraped,1)-2;
sfigure(4882);
clf
errorbar(bc_angles_wraped(1:hilight_cut),to_vals_lin_quad(1:hilight_cut,2)-lin_fit_max(1),to_vals_lin_quad(1:hilight_cut,3),'ko')
hold on
errorbar(bc_angles_wraped(hilight_cut+1:end),to_vals_lin_quad(hilight_cut+1:end,2)-lin_fit_max(1),to_vals_lin_quad(hilight_cut+1:end,3),'bx')
xlabel('pol angle')
ylabel('to value')
set(gcf,'color','w')
%set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1600, 900])



%% Plot the purity data
hilight_cut=size(bc_angles_wraped,1)-2;
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

sfigure(6219);
clf
errorbar(hwp_ang(is_outlier_idx),vec_corr_to(is_outlier_idx)-lin_fit_max(1),to_vals_lin_quad(is_outlier_idx,3),'ko')
hold on
errorbar(hwp_ang(~is_outlier_idx),vec_corr_to(~is_outlier_idx)-lin_fit_max(1),to_vals_lin_quad(~is_outlier_idx,3),'rx')
xlabel('hwp angle')
ylabel('to value')
set(gcf,'color','w')
%set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1600, 900])

li_yan = [413.085836,413.084856,413.083876,413.082896];
li_yan_pol = [0, 35.26, 54.74, 90];
hebec_constants
li_yan_freq = const.c./li_yan;