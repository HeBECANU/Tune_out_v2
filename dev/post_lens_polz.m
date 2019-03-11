
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

naked_pol1 = [0:20:180;
            364,205,220,233,300,305,220,201,212,300;
            40,18,196,297,243,236,206,2,42,33;
            6.45,6.4,57.6,133,59,0.3,40,87,47.5,4.5;
            138,118,91,148,324,307,289,104,152,135]';
        
naked_pol2 = [10:20:170;
                336,230,160,173,270,248,286,250,300;
                227,196,183,47,237,237,2,259,41;
                0.25,21,82.5,101.7,13,10,99,96.5,39;
                126,98,85,340,315,296,280,146,322]';
naked_pol3 = [55,65,135,145,175;
            251,258,207,219,326;
            332,271,202,258,231;
            154,100,158,194,23;
            nan,158,261,318,318]';
naked_pol = [naked_pol1;naked_pol2;naked_pol3];
[~,idx] = sort(naked_pol(:,1));
naked_pol = naked_pol(idx,:);

d_th_max = naked_pol(:,1) - naked_pol(:,3);
d_th_min = naked_pol(:,1) - naked_pol(:,5);
naked_power = naked_pol(:,2)+naked_pol(:,4);

contrast = @(pol_data) (pol_data(:,2) - pol_data(:,4))./(pol_data(:,2) + pol_data(:,4));
ratio = @(pol_data) pol_data(:,4)./pol_data(:,2);
q = @(p) p(:,1);
qwpa = 0:20:18;

sfigure(5831);
clf
subplot(2,2,1)
plot(q(pol_pre_comp),contrast(pol_pre_comp),'kx-')
hold on
plot(q(pol_post_comp),contrast(pol_post_comp),'rx-')
plot(q(pol_post_aptr),contrast(pol_post_aptr),'ro-')
plot(q(naked_pol),contrast(naked_pol),'go-')
ylim([0,1])
xlabel('QWP angle')
ylabel('Contrast')
legend('Pre-window','After chamber','After chamber with aperture','At window','location','SouthWest')
title('Contrast pre and post window')
subplot(2,2,2)
plot(q(pol_pre_comp),ratio(pol_pre_comp),'kx-')
hold on
plot(q(pol_post_comp),ratio(pol_post_comp),'rx-')
plot(q(pol_post_aptr),ratio(pol_post_aptr),'ro-')
plot(q(naked_pol),ratio(naked_pol),'go-')
ylim([0,1])
xlabel('QWP angle')
ylabel('min/max ratio')
title('Power ratio pre and post window')

subplot(2,2,3)
plot(q(naked_pol),naked_power,'gx-')
ylim([0,max(naked_power)])
title(sprintf('Total power, [mean,std] = [%.2f,%.2f]',mean(naked_power),std(naked_power)))
xlabel('QWP angle')
ylabel('Pmax+Pmin')


subplot(2,2,4)
plot(q(naked_pol),d_th_min,'b*-')
hold on
plot(q(naked_pol),d_th_max,'r*-')
xlabel('QWP angle')
ylabel('Polarizer angle offset')
legend('Minimizing angle','Maximizing angle','location','NorthWest')