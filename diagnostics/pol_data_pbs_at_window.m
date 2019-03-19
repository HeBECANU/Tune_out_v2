fwtext('Begin')
% 
% background=0.1;%allpowerinmicrowatts
%Anglereadingsondialwhenfastaxisvertical
% NB test HWP mounted with notch horizontal
hwp_offset=20;
qwp_offset=24;
% Calibrating fast axes:
% hwp 337 qwp 198 min transmission of 1.8 to reflection of 490
% with test hwp in ^ position 155 min 0.99 to 490: VERT IN, FAST AXIS VERT?
%QWP HWP QWP_anal HWP_Anal :155 Phi_min Phi_max (Transmit,Reflect)...
data=  [130,333,nan,155,343,140,332,131,112,108,347,196;
        134,333,nan,155,300,162,325,160,242,130,352,286;
        138,333,nan,155,273,185,275,180,068,150,314,109;
        142,333,nan,155,252,202,256,202,169,178,285,209;
        146,333,nan,155,245,216,270,185,270,156,307,124;
        150,333,nan,155,242,222,297,165,297,138,335,130;
        154,333,nan,155,235,256,325,138,079,116,380,117;
        162,333,nan,155,195,278,380,66.,039,084,409,084;
        177,333,nan,155,087,282,527,10.2,122,12.,421,165;
        187,333,nan,155,022,406,500,03.6,204,0.8,437,250;
        202,333,nan,155,21.7,412,420,38.8,287,27.8,418,242;
        220,333,nan,155,150,306,350,129,275,125,367,232;
        226,333,nan,155,199,262,309,157,179,159,352,126;
        234,333,nan,155,240.,210,260,208,172,220,282,208;
        246,333,nan,155,242.,211,353,158,226,147,346,272;
        254,333,nan,155,267.,222,411,111,037,100,395,351;
        280,333,nan,155,021.,477,535,03.7,204,01.3,508,249;
        283,333,nan,155,13.7,464,540,05.2,204,04.6,463,158;
        270,340,nan,155,126.,360,412,105,121,076,368,77;
        274,350,nan,155,203.,280,317,207,024,184,291,86;
        268,354,nan,155,299.,153,297,155,066,160,312,21;
        280,010,nan,155,472,045,465,41.7,066,34,437,22;
        290,006,nan,155,326,158,372,148,059,140,366,14];
%         270,333,nan,155,080.,415,454,24.3,299,15.2,489];

% QWP rotation data
% Set transmission to ZERO with control waveplates:
% HWP 283 fast horizontal
% QWP 333 fast 45deg to horizontal
% PBS powers 287,279
% Insert test HWP and rotate to 302 degrees, 
% PBS powers 17, 517
% contrl QWP, HWP, power 1, power 2
qwp_data = [130,333,405,95
        134,333,460,76
        138,333,477,53
        142,333,410,35
        146,333,417,29.6
        150,333,437,32.6
        154,333,408,36
        162,333,374,68
        177,333,292,192
        187,333,191,297
        202,333,76,420
        220,333,15.7,482
        226,333,14,464
        234,333,23,460
        246,333,42,440
        254,333,62,421
        270,333,116,370
        280,333,175,333
        283,333,202,302
        270,340,48,440
        274,350,20.5,470
        268,354,64,445
        280,010,211,313
        290,006,77,428
        198,337,157,345];

    
    
sfigure(412);
clf
qm = qwp_data(:,2) == 333;
subplot(2,1,1)
plot(qwp_data(qm,1), qwp_data(qm,4))
hold on
plot(qwp_data(qm,1), qwp_data(qm,3))
subplot(2,1,2)
plot(qwp_data(qm,1), contrast(qwp_data(qm,3:4)),'*')
% title('Contrast 
    
    
QWP = data(:,1);
HWP = data(:,2);
QWT = data(:,3);
th_ref = data(:,4);
p_ref = data(:,5:6);
P_1 = data(:,7:8);
th_1 = data(:,9);
P_2 = data(:,10:11);
th_2 = data(:,12);

HWmask = HWP == 333;
% Max observed contrast is the mod-square of the projection onto the Poincare plane spanned by
% linear polz bases, so circular amplitude is sqrt(1-contrast), with sign left undetermined
contrast = @(x) abs(diff(x,1,2)./sum(x,2));
ratio = @(x) min(x,[],2)./max(x,[],2);
circmdl = @(p,x) abs(p(1)*cos(p(3)*x/180+p(2)))+p(4);
beta0 = [1,deg2rad(140),5,0];
fit1 = fitnlm(QWP_corr(HWmask),sqrt(1-contrast(P_1(HWmask,:)).^2),circmdl,beta0);
fit2 = fitnlm(QWP_corr(HWmask),sqrt(1-contrast(P_2(HWmask,:)).^2),circmdl,beta0);


th_test = mod(2.*[th_ref,th_1,th_2],180);
th_test = deg2rad(th_test);
%mod(pi/2-2.*th_test(:,2),2*pi)
p_contrast = abs([p_ref(:,1),P_1(:,1),P_2(:,1)]-[p_ref(:,2),P_1(:,2),P_2(:,2)])./([p_ref(:,1),P_1(:,1),P_2(:,1)]+[p_ref(:,2),P_1(:,2),P_2(:,2)]);
contrast_mdl = @(b,x) abs(b(1).*sin(2.*x(:,1)+b(2)));
beta = [0.1,3.3;
        0.4,1.0;
        0.3,0.3;
        0.2,2.2;
        0.1,0.3;
        0.1,1.4;
        0.8,1.9;
        0.7,4.0;
        0.1,5.0;
        0.1,2.0;
        0.1,5.5;
        0.5,3.5];
    
    %%
stokes_data = zeros(size(th_test,1),2);
stokes_unc = zeros(size(th_test,1),2);
t = linspace(0,2.*pi);
for ii = 1:size(th_test,1)
    unc_tol = 1e9;
    beta(8,:) = [(contrast(P_1(ii,:))+contrast(P_2(ii,:)))./2,mod(pi/2-2.*th_test(ii,1),2*pi)];
    for jj = 1:size(beta,1)
        fit_ratios = fitnlm(th_test(ii,:)',p_contrast(ii,:)',contrast_mdl,beta(jj,:));
        if fit_ratios.Coefficients.SE(1)<unc_tol
            unc_tol = fit_ratios.Coefficients.SE(1);
            fit_vals = fit_ratios.Coefficients.Estimate;
            fit_unc = fit_ratios.Coefficients.SE;
        end
    end
    stokes_data(ii,:) = fit_vals;
    stokes_unc(ii,:) = fit_unc;
    sfigure(55995959);
    clf
    scatter(th_test(ii,:)',p_contrast(ii,:)','*')
    hold on
    plot(t,contrast_mdl(fit_ratios.Coefficients.Estimate,t'))
    title(['Control QWP =',num2str(QWP(ii))])
    xlabel('Test HWP angle')
    ylabel('contrast')
    pause(3.5)
end



sfigure(105289375);
clf

QWP_corr = mod(QWP,360);
Theta = linspace(min(QWP_corr),max(QWP_corr),150);
subplot(6,1,[1,2])
plot(QWP_corr(HWmask),sqrt(1-contrast(P_1(HWmask,:)).^2),'*')
hold on
plot(QWP_corr(HWmask),sqrt(1-contrast(P_2(HWmask,:)).^2),'*')
plot(Theta,circmdl(fit1.Coefficients.Estimate,Theta))
plot(Theta,circmdl(fit2.Coefficients.Estimate,Theta))
scatter(QWP_corr(HWmask),sqrt(1-stokes_data(HWmask,1).^2),'*')
% plot(QWP_corr(~HWmask),contrast(P_2(~HWmask,:)),'*')
% plot(QWP_corr(HWmask),contrast(P_2(HWmask,:)),'*')
xlabel('Control QWP angle')
ylabel('Circular at measured angles')
legend('Measurement 1','Measurement 2','Fit 1','Fit 2','Multi Fit','Rotating control HWP','location','SouthWest')
ylim([0,1])
title('|Circular polarization component|?')

% % Hm, not super sure about the interpretation of subplot 1. Should double check the expressions
% that come out. Cf subplot 2, which shows  nice |sin(x)| behaviour (but slightly less than perfect
% contrast, might have been mis-re-alignment of HWP angle...
% 

subplot(6,1,[3,4])

plot(QWP_corr(HWmask),1-contrast(p_ref(HWmask,:)),'*')
hold on
xlabel('Control QWP angle')
ylabel('Circular component')
title('Contrast at reference angle')


subplot(6,1,5)
plot(QWP_corr(HWmask),P_1(HWmask,:),'*')
hold on
plot(QWP_corr(HWmask),P_2(HWmask,:),'*')
title('Raw power measurements')

subplot(6,1,6)
plot([120,300],[45,45],'k-.')
hold on
plot(QWP,mod(th_1-th_2,180)-90,'r*')
plot([120,300],[-45,-45],'k-.')
title('Angle difference between measurements')
ylabel('\Delta \theta_{HWP}-\pi/2')
xlabel('QWP angle')



%% Some test from K
mdl = @(p,x) abs(p(1).*cos(p(3).*x.*pi/180+p(2)));
beta0 = [1.0,0.5,2.0];
wlin = 1./stokes_unc(HWmask,1).^1;
wlin = wlin./sum(wlin);
contrast_fit = fitnlm(QWP_corr(HWmask),stokes_data(HWmask,1),mdl,beta0,'Weights',wlin);
contrast_p = (contrast(P_1(HWmask,:))+contrast(P_2(HWmask,:)))./2;
% contrast_fit = fitnlm(QWP_corr(HWmask),contrast_p,mdl,beta0);

sfigure(1231231);
clf
ci_size_disp = 1-erf(1/sqrt(2));
x_samp_pad = linspace(0,360,1e5);
[y_samp_val,y_samp_ci]=predict(contrast_fit,x_samp_pad','Alpha',ci_size_disp); %
patch([x_samp_pad, fliplr(x_samp_pad)], ([y_samp_ci(:,1)', fliplr(y_samp_ci(:,2)')]), color_shaded,'EdgeColor','none');  %
hold on
scatter(QWP_corr(HWmask),contrast_p,'*')
plot(x_samp_pad,y_samp_val,'-','color',colors_main(3,:),'LineWidth',1.5)
errorbar(QWP_corr(HWmask),stokes_data(HWmask,1),stokes_unc(HWmask,1),'kx')
ylim([0,1])
xlabel('QWP angle')
ylabel('Contrast')
%%
sfigure(459640)
clf
V_run = pol_data_post(:,7).*2.*sqrt(pol_data_post(:,3).*pol_data_post(:,5))...
./(pol_data_post(:,3)+pol_data_post(:,5));%the V parameter for each run
scatter(pol_data_post(38:end-5,1),V_run(38:end-5),'*')
hold on
A = sqrt(1-stokes_data(HWmask,1).^2).*(1-2.*(abs(QWP_corr(HWmask)-234)<45));
scatter(QWP_corr(HWmask),A,'kx')
A_2 = sqrt(1-contrast_p.^2).*(1-2.*(abs(QWP_corr(HWmask)-234)<45));
scatter(QWP_corr(HWmask),A_2,'ro')
xlabel('QWP angle')
ylabel('Circular component A')
legend('Post','Pre','location','Best')
title('Comparison between Pre and Post')



fwtext('Done!')

%% checking polz stability; in short,looks pretty good
% dir='C:\Users\BEC Machine\Dropbox\MATLAB\Tune_out_v2_trap_freq_master\diagnostics\pol_stab';
% f1='DATA27.csv';
% f2='DATA28.csv';
% in1=importpowermeter(fullfile(dir,f1));
% in2=importpowermeter(fullfile(dir,f2));
% 
% T1=in1.time;
% P1=in1.data;
% T2=in2.time;
% P2=in2.data;
% M1=P1/max(P1)>1e-1;
% M2=P2/max(P2)>1e-2;
% Tf1=T1(M1);
% Tf2=T2(M2);
% Pf1=P1(M1);
% Pf2=P2(M2);
% 
% P_sum = Pf1+Pf2;
% contrast=(Pf2-Pf1)./P_sum;
% C_err = contrast-mean(contrast);
% 
% sfigure(784);
% clf;
% subplot(4,2,[1,2])
% plot(Tf1,Pf1)
% hold on
% plot(Tf2,Pf2)
% plot(Tf2,P_sum)
% %ylim(.06,.09)
% xlabel('Time(s)')
% ylabel('Power (W)')
% title('Power at each meter')
% legend('Transmitted','Reflected','total')
% 
% subplot(4,2,[3,4])
% plot(Tf1,Pf1/max(Pf1))
% hold on
% plot(Tf2,Pf2/max(Pf2))
% title('Scaled power at each meter')
% legend('Transmitted','Reflected')
% xlabel('time (ms)')
% ylabel('Power/max(P) ')
% 
% subplot(4,2,[5,6])
% plot(contrast)
% % ylim([0,1])
% title('Contrast vs time')
% xlabel('time (ms)')
% ylabel('contrast')
% 
% subplot(4,2,[7,8])
% histogram(C_err,100,'Normalization','PDF');
% title('Contrast fluctuations')
% xlabel('Deviation from mean')
% ylabel('Probability')
% % 
% % subplot(4,2,8)
% % plot(Pf1,Pf2,'*')
% % % plot(Pf1-mean(Pf1),C_err,'*')
% % hold on
% % % plot(Pf2-mean(Pf2),C_err,'*')
% % xlabel('Total power variation')
% % ylabel('Contrast variation')
% % title('Power and contrast variations')
% % legend('transmitted power','reflected power')