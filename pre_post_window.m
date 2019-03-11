
% % Define data
%pre window
pol_data_pre = [280,45.6,171,1.2,265,46;
            310,34.2,153,8.5,247,46;
            270,39.9,153,8.5,265,46;
            220,42.2,40,1.6,305,45;
            150,43.2,139,2.24,47,47;
            130,36.9,159,7.5,247,45
            ];

pre_pows = pol_data_pre(:,2) + pol_data_pre(:,4);
pre_frac = pol_data_pre(:,2)./pre_pows;
%post window        
pol_data_post = [280,168,257,0.13,171,nan;
            310,140.7,117,40.8,205,nan;
            270,167.3,248,5.37,159,nan;
            220,110,273,44.1,192,nan;
            150,90.6,245,63.9,152,nan;
            130,119,284,31,12,nan
            ];
post_pows=pol_data_post(:,2) + pol_data_post(:,4);
post_frac = pol_data_post(:,2)./post_pows;


pol_data_pre_long = [280,213,193,3.3,102.5,230;
                    310,151,38,54,117,210;
                    270,209,15,32.6,101,269;
                    220,233,152.5,11.8,241,250;
                    150,247,221,3.1,133,268;
                    130,186,35.5,54.5,119,254;
                    142,175,46,34.2,140,245;
                    162,195,50,0.69,322,204;
                    202,210,176,9.2,84,220;
                    226,178,146,12.7,246,230;
                    254,120,342,82,67,233;
                    210,212,163,8.47,72,215];
pol_dat_quick = [0,30,60,90,120,150,180,210,240,270,300,330;
                280,266,213,233,239,250,256,236,205,223,235,263;
                10.6,14.1,57.6,33.3,26.8,21,10.2,9.7,53,41.2,35,20.9]';
pol_quick_plus = [15,45,75,105,135,165,195,225,255,285,315,345;
                199,239,197,271,221,290,266,245,163,280,220,290;
                10.5,25,91.8,0.6,62,0.6,13.3,17,116,0.75,64,57]';
pol_dat_late = [pol_dat_quick;pol_quick_plus];


pol_pre_comp = [0,30,45,60;
                190,230,240,197;
                30,338,330,324;
                7.9,10.8,25,46;
                122,78,240,236]';
pol_post_comp = [0,30,45,60;
                290,230,206,163;
                60,269,92,224;
                10.5,40,117,100;
                335,181,24,304]';
    

pol_post_aptr = [0,30,45,60;
                122,139,94,97;
                246,270,114,220;
                4.6,21.4,57,55;
                156,3,14,121]';

PBS_pol = [0:20:340;
            246,270,253,195,194,237,243,235,274,253,247,225,209,189,293,256,249,307;
            10.7,9.4,28.7,72,85,7.9,32.8,59,4.5,13.2,16.2,18,68,117,9.3,38.9,66,4.5]';


% Quick functions
contrast = @(pol_data) (pol_data(:,2) - pol_data(:,4))./(pol_data(:,2) + pol_data(:,4));
ratio = @(pol_data) pol_data(:,4)./pol_data(:,2);
quick_contrast = @(pol_data) (pol_data(:,2) - pol_data(:,3))./(pol_data(:,2) + pol_data(:,3));


%%
sfigure(42);
clf;
subplot(2,3,1)
A_post = 2.*sqrt(pol_data_post(:,2).*pol_data_post(:,4))./(pol_data_post(:,4)+pol_data_post(:,2));
sign_flip = ones(size(pol_data_post,1),1).*(1-2*(45>abs(pol_data_post(:,1)-235)));
A_post = A_post.*sign_flip;
A_pre = 2.*sqrt(pol_data_pre(:,2).*pol_data_pre(:,4))./(pol_data_pre(:,4)+pol_data_pre(:,2));
A_pre_long = 2.*sqrt(pol_data_pre_long(:,2).*pol_data_pre_long(:,4))./(pol_data_pre_long(:,4)+pol_data_pre_long(:,2));

scatter(pol_data_pre(:,1),A_pre,30,'kx')
hold on
scatter(pol_data_post(:,1),A_post,30,'bo')
scatter(pol_data_pre_long(:,1),A_pre_long,30,'r*')
xlabel('QWP')
ylabel('A')
legend('pre','post','pre long','location','southEast')
set(gcf,'color','w')
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1600, 900])

subplot(2,3,2)
scatter(pol_data_pre(:,1),pol_data_pre(:,4)./(pol_data_pre(:,2)),30,'kx')
hold on
scatter(pol_data_post(:,1),pol_data_post(:,4)./(pol_data_post(:,2)),30,'bo')
scatter(pol_data_pre_long(:,1),pol_data_pre_long(:,4)./(pol_data_pre_long(:,2)),30,'r*')
xlabel('QWP')
ylabel('p_{min}/p_{max}')
legend('pre','post','pre long','location','southEast')
set(gcf,'color','w')
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 1600, 900])

subplot(2,3,3)
scatter(pol_data_pre(:,1),pol_data_pre(:,4)./(pol_data_pre(:,4)+pol_data_pre(:,2)))
hold on
scatter(pol_data_post(:,1),pol_data_post(:,4)./(pol_data_post(:,4)+pol_data_post(:,2)))
hold off
xlabel('QWP')
ylabel('p_{min}')
legend('pre','post','location','southEast')

subplot(2,3,4)
scatter(pol_data_pre(:,1),pol_data_pre(:,2)./(pol_data_pre(:,4)+pol_data_pre(:,2)))
hold on
scatter(pol_data_post(:,1),pol_data_post(:,2)./(pol_data_post(:,4)+pol_data_post(:,2)))
hold off
xlabel('QWP')
ylabel('p_{min}')
legend('pre','post','location','southEast')

subplot(2,3,5)
scatter(pol_data_pre(:,1),2.*sqrt(pol_data_pre(:,2).*pol_data_pre(:,4))./(pol_data_pre(:,4)+pol_data_pre(:,2)))
hold on
scatter(pol_data_post(:,1),2.*sqrt(pol_data_post(:,2).*pol_data_post(:,4))./(pol_data_post(:,4)+pol_data_post(:,2)))
hold off
xlabel('QWP')
ylabel('A')
legend('pre','post','location','southEast')

subplot(2,3,6)
plot(pol_data_pre(:,1),pre_frac,'kx')
hold on
plot(pol_data_post(:,1),post_frac,'bo')
% errorbar(pol_data_pre_long(:,1),pre_frac_long,pow_errs,'r*')
legend('Pre','Post','Janky table')

sfigure(44);
clf;
subplot(2,1,1)
plot(pol_dat_late(idx,1),quick_contrast(pol_dat_late(idx,:)),'x-')
hold on
plot(PBS_pol(:,1),quick_contrast(PBS_pol),'r-.')
plot(pol_data_pre(:,1),contrast(pol_data_pre),'kx')
plot(pol_data_pre_long(:,1),contrast(pol_data_pre_long),'rx')
plot(pol_data_post(:,1),contrast(pol_data_post),'bo')

xlabel('QWP angle')
ylabel('(P_{max}-P_{min})/P_{total}')
legend('Janky table late night','Janky table PB cube','Before window','Janky table','After chamber','location','SouthWest')
title('Contrast')
subplot(2,1,2)
plot(pol_data_pre(:,1),mod(pol_data_pre(:,3)-pol_data_pre(:,5)+45,90)-45,'kx')
hold on
plot(pol_data_pre_long(:,1),mod(pol_data_pre_long(:,3)-pol_data_pre_long(:,5)+45,90)-45,'rx')
plot(pol_data_post(:,1),mod(pol_data_post(:,3)-pol_data_post(:,5)+45,90)-45,'bo')
xlabel('QWP angle')
ylabel('mod(\theta_{max}-\theta_{min},\pi/2)')
title('Polarizer angle error from 90\circ')
legend('Before window','Janky table','After chamber')

sfigure(686);
[~,idx] = sort(pol_dat_late(:,1));
plot(pol_dat_late(idx,1),quick_contrast(pol_dat_late(idx,:)),'x-')
xlabel('QWP angle')
ylabel('(Pmax-Pmin)/(Pmax+Pmin)')
title('Contrast versus QWP on jank table')


% % Comparing pre & post including aperture
sfigure(5830);
subplot(1,2,1)
plot(qwpa,contrast(pol_pre_comp),'kx-')
hold on
plot(qwpa,contrast(pol_post_comp),'rx-')
plot(qwpa,contrast(pol_post_aptr),'ro-')
ylim([0,1])
xlabel('QWP angle')
ylabel('Contrast')
legend('Pre-window','After chamber','After chamber with aperture','location','SouthWest')
title('Contrast pre and post window')
subplot(1,2,2)
plot(qwpa,ratio(pol_pre_comp),'kx-')
hold on
plot(qwpa,ratio(pol_post_comp),'rx-')
plot(qwpa,ratio(pol_post_aptr),'ro-')
ylim([0,1])
xlabel('QWP angle')
ylabel('min/max ratio')
title('Power ratio pre and post window')


% % plot all the things
sfigure(3333);
plot(pol_dat_late(idx,1),quick_contrast(pol_dat_late(idx,:)),'x-')
hold on
plot(PBS_pol(:,1),quick_contrast(PBS_pol),'r-.')
plot(pol_data_pre(:,1),contrast(pol_data_pre),'kx')
plot(pol_data_pre_long(:,1),contrast(pol_data_pre_long),'rx')
plot(pol_data_post(:,1),contrast(pol_data_post),'bo')
plot(qwpa,contrast(pol_pre_comp),'kx-')
plot(qwpa,contrast(pol_post_comp),'gx-')
plot(qwpa,contrast(pol_post_aptr),'go-')
%% Attempting to recover post from pre
del = 0.15705;%0.22034;%0.19263;
phi = 0.22308;%-0.203;
A=2.*sqrt(pol_data_pre(:,2).*pol_data_pre(:,4))./(pol_data_pre(:,4)+pol_data_pre(:,2));
psi = pol_data_pre(:,5);
S_3 = A.*[-1;1;-1;-1;-1;1];
S_2 = cos(2.*psi).*cos(asin(A));
S_1 = sin(2.*psi).*cos(asin(A));

A_post = 2.*sqrt(pol_data_post(:,2).*pol_data_post(:,4))./(pol_data_post(:,4)+pol_data_post(:,2));
A_post = A_post.*[-1;1;-1;-1;-1;1];
A_recon = sin(del).*(S_2.*cos(psi.*pi./90+phi.*2*pi)-S_1.*sin(psi.*pi./90+phi.*2*pi))+S_3.*cos(del);
mdl = @(b,x) sin(b(1)).*(x(:,2).*cos(x(:,4).*pi./90+b(2).*2*pi)-x(:,1).*sin(x(:,4).*pi./90+b(2).*2*pi))+x(:,3).*cos(b(1));
beta0 = [0.19263,0.0];
fit_mdl = fitnlm([S_1,S_2,S_3,psi],A_post,mdl,beta0);
sfigure(47);
clf
scatter(pol_data_pre(:,1),A_recon)
hold on
scatter(pol_data_post(:,1),A_post)
%legend('reconstructed A','measured post A')
hold off
xlabel('QWP')
ylabel('A')

