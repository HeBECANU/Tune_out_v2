phi = 0.4;%0.628547324144821;%3.1227;%
%qwp data
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
     270,725733518.4,83,725733480.3,96;
     246,725729215.0,103,725729418,121;
     283,725735078.9,45,725735160.3,81;
     283,725736295.0,22,725736324.0,35;
     280,725735737.3,36,725735760.8,54
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
   270,167.3,248,5.37,159;
   246,127,29,45,136;
   283,195,83,0.6,352;
   283,195,83,0.6,352;
   280,168,257,0.13,171
   ];%
sign_flip = ones(size(polz_data ,1),1).*(1-2*(45>abs(polz_data(:,1)-235)));
A=2*sqrt(polz_data(:,2).*polz_data(:,4))./(polz_data(:,2)+polz_data(:,4)).*sign_flip;
p_min = polz_data(:,4)./(polz_data(:,4)+polz_data(:,2));
p_max = polz_data(:,2)./(polz_data(:,4)+polz_data(:,2));
P=p_min-(p_min-p_max).*cos(polz_data(:,5).*pi/180+2.*pi.*phi).^2;
TO = to_vals_lin_quad(:,2);
psi = polz_data(:,5);
qwp_num=size(to_vals_lin_quad,1);

%hwp data
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
    0.34,111.0,79.0,20.0,NaN,NaN,304,1;
    0.37,130.5,120.0,40.0,nan,nan,314,1;
    0.14,151.0,98.0,63.0,nan,nan,325.0,1;
    0.13,170,113,84,nan,nan,335,0;
    0.68,190,133,104,nan,nan,344,0;
    0.94,29,89,120,85,119,354,0;
    0.96,50.5,75,146,nan,nan,6,0;
    0.59,70,120,165,nan,nan,15,0;
    0.06,92,77,181,nan,nan,26,1;
    0.52,61,84,154,nan,nan,10,0;
    0.86,46,74,137,nan,nan,3,0;
    0.54,24.5,72,112,nan,nan,352,0;
    0.19,180,93,89,nan,nan,339,0;
    0.12,160,120,70.5,nan,nan,319,1;
    0.33,140.5,112,50,118,50,320,1;
    0.36,121,86,32,nan,nan,310,1;
    0.17,100,72,191,nan,nan,300,1;
    0.07,80,107.5,165,nan,nan,290,0;
    0.25,99.5,72,188,72,187.5,299,1;
    0.37,120,78,221,nan,nan,309,1;
    0.24,140,62,230,nan,nan,320,1;
    0.23,147,87.6,238,nan,nan,324,1;
    0.11,155,60.4,248,nan,nan,327,1;
    0.14,160,76,249,nan,nan,329,1;
    0.20,165,80,255,77,254,332,1;
    0.16,171,79,264,nan,nan,335,0;
    0.23,176.5,65,266,nan,nan,338.5,0;
    0.23,176.5,65,266,nan,nan,338.5,0;
    0.33,181,71,271,nan,nan,340,0;
    0.52,187,92,278,nan,nan,342,0;
    0.66,194,78,282,nan,nan,347,0;
    0.72,199,82,290,80,291,350,0;
    0.65,208,57,298,nan,nan,353,0;
    0.90,217,78,308,nan,nan,358,0;
    0.70,231,70.1,135,nan,nan,4,0;
    0.64,239.5,53,150,60,151,9.5,0;
    0.58,249,75.1,160,nan,nan,14.5,0];
polz_sign=polz_data(:,8);
polz_sign(~polz_sign)=-1;
A = [A;2*sqrt(polz_data(:,1).*polz_data(:,3))./(polz_data(:,1)+polz_data(:,3)).*-polz_sign];
p_min = polz_data(:,1)./(polz_data(:,1)+polz_data(:,3));
p_max = polz_data(:,3)./(polz_data(:,1)+polz_data(:,3));
P=[P;p_min-(p_min-p_max).*cos(polz_data(:,2).*pi/180+2.*pi.*phi).^2];
TO = [TO;to_vals_lin_quad(:,2)];
psi = [psi;polz_data(:,2)];

%%
polz_sign=polz_data(:,8);
polz_sign(~polz_sign)=-1;
vec_amp=-6626.8;%0;%-6e3;%19e3;
A_temp = 2*sqrt(polz_data(:,1).*polz_data(:,3))./(polz_data(:,1)+polz_data(:,3));
vec_corr_to=to_vals_lin_quad(:,2)-A_temp.*polz_sign*vec_amp;
modelfun = @(b,x) b(1).*(x(:,3)+x(:,2).*cos(x(:,1).*pi./180+b(2).*2*pi).^2)+b(3);
opts = statset('MaxIter',1e5);
beta0 = [1e3,0.6,nanmean(vec_corr_to)]; %intial guesses
wlin=1./(to_vals_lin_quad(:,3).^2);
fit_mdl_lin = fitnlm([polz_data(:,2),p_max-p_min,p_min],vec_corr_to,modelfun,beta0,...
    'Options',opts,'Weights',wlin,'CoefficientNames' ,{'amp','phase','offset'});
lin_fit_max=[fit_mdl_lin.Coefficients.Estimate(1)+fit_mdl_lin.Coefficients.Estimate(3)...
    sqrt(fit_mdl_lin.Coefficients.SE(1)^2+fit_mdl_lin.Coefficients.SE(3)^2)];

sfigure(2171);
clf
ci_size_disp=1-erf(1/sqrt(2));%one sd %confidence interval to display
plot_padd=20;
xsamp = linspace(min(polz_data(:,2))-plot_padd,max(polz_data(:,2))+plot_padd,1e4).';
fit_values = fit_mdl_lin.Coefficients{:,1};
hold on
title(['amp= ',num2str(fit_values(1)),'MHz , phase= ',num2str(fit_values(2)),', offset= ',num2str(fit_values(3)),'MHz , period=180(fixed)^\circ'])
xlabel('Input Pol angle (degrees)')
ylabel(sprintf('Tune out value-%.1f±%.1f (MHz)',lin_fit_max(1),lin_fit_max(2)) )
fprintf('Tune out value-%.1f±%.1f (MHz) (BLUE)\n',lin_fit_max(1),lin_fit_max(2))
fprintf('Tune out value-%.1f±%.1f (MHz) (RED)\n',lin_fit_max(1)/2,lin_fit_max(2)/2)
errorbar(polz_data(:,2),vec_corr_to-lin_fit_max(1),to_vals_lin_quad(:,3),'ko')
% [y_lin,yci_lin]=predict(fit_mdl_lin,xsamp,'Prediction','observation','Alpha',ci_size_disp); %'Prediction','observation'
y_lin =fit_values(1).*(cos(xsamp.*pi/180+fit_values(2).*2*pi).^2)+fit_values(3);
plot(xsamp,y_lin-lin_fit_max(1),'b-','LineWidth',1.6)
% plot(xsamp,yci_lin(:,1)-lin_fit_max(1),'r-','LineWidth',1.6)
% plot(xsamp,yci_lin(:,2)-lin_fit_max(1),'r-','LineWidth',1.6)

set(gcf,'color','w')

%A = A.*(1-2.*(TO<nanmean(TO)));
%%
mdl = @(b,x) b(1).*x(:,1)+b(2).*x(:,2)+b(3);
beta0 = [2000,100,725735959.15];
fit_mdl = fitnlm([A,P],TO,mdl,beta0);

A_samp = linspace(-1,1)';
P_samp = linspace(0,1)';
[A_samp,P_samp] = meshgrid(A_samp,P_samp);
b = fit_mdl.Coefficients{:,1};
T_samp = b(1).*A_samp+b(2).*P_samp+b(3);

TO_res = TO-(b(1).*A+b(2).*P+b(3));

sfigure(4008);
clf
offset = nanmean(TO);
scatter3(A,P,TO-offset)
hold on
s = surface(A_samp,P_samp,T_samp-offset);
alpha(s,'z')
xlabel('A')
ylabel('P')
sfigure(4009);
scatter3(A,P,TO_res)
xlabel('A')
ylabel('P')
%%
mean(abs(TO_res))
fit_mdl
%%
figure(11109);
[x y z] = sphere(256);
h = surfl(x, y, z);
colormap(viridis)
set(h, 'FaceAlpha', 0.7)
shading interp
hold on
S_3 = A;
S_2 = cos(2.*psi).*cos(asin(A));
S_1 = sin(2.*psi).*cos(asin(A));
scatter3(S_1(1:qwp_num),S_2(1:qwp_num),S_3(1:qwp_num),'bx')
scatter3(S_1(qwp_num+1:end),S_2(qwp_num+1:end),S_3(qwp_num+1:end),'r*')
hold off
axis equal
title('Poincare Sphere Representation')
xlabel('S_1')
ylabel('S_2')
zlabel('S_3')