fwtext('Running polz analysis')

% Fun defns


% Data input

offset = -.0;
qwp = [280,283,246,270,286,310,254,234,226,220,202,187,177,162,154,130,134,138,142,146,150,270,274,280,290];
hwp = [333,333,333,333,333,333,333,333,333,333,333,333,333,333,333,333,333,333,333,333,333,340,350,10,6];
pmax = [230,233,259,245,300,238,170,219,111,205,183,255,239,220,199,165,147,143,191,175,177,216,199,200,205];
thmax = [217,217,255,234,21,4,40,344,188,7,209,27,241,242,252,294,198,210,5,25,39,212,20,122,295];
pmin = [0.2,0.8,60,6,2.1,54,55,77,56,34,6.6,.77,16.5,30.5,110,90,72.1,76,98,94.9,126,57,81.5,22,120];
thmin = [305,123,140,139,301,100,142,116,277,107,116,308,319,151,332,271,99,104,129,141,330,312,125,214,30];
pdata = [qwp;hwp;pmax;thmax;pmin;thmin]';


post_data = [
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
   310,140.7,117,40.8,205;
   286.5,165,88.5,2.18,179;
   270,167.3,248,5.37,159;
   246,127,29,45,136;
   283,195,83,0.6,352;
   283,195,83,0.6,352;
   280,168,257,0.13,171
   ];%

post_hwp = 333*ones(size(post_data,1),1);
post_data = [post_data(:,1),post_hwp,post_data(:,2:end)];

% Compute things

% Plotting
bigplot(pdata,offset,'Pre window',1)
bigplot(post_data,offset,'Post window, HWP 333',2)




corrplot(pdata,offset,'Correlations pre window',3)
corrplot(post_data,offset,'Correlations post window',4)

d_th_max = wrapToPi(deg2rad(pdata(:,1) - pdata(:,3)));
d_th_min = wrapToPi(deg2rad(pdata(:,1) - pdata(:,5)));
post_power = pow_sum(post_data);
pre_power = pow_sum(pdata);

[~,pre_idx] = sort(pdata(:,1));
[~,post_idx] = sort(post_data(:,1));

pdata = pdata(pre_idx,:);
post_data = post_data(post_idx,:);

sfigure(578);
clf;
subplot(3,1,1)
plot(q(pdata),contrast(offset,pdata),'k*-')
hold on
plot(q(post_data),contrast(offset,post_data),'r*-')
ylim([0,1])
xlabel('QWP angle')
ylabel('Contrast')
title('Contrast')

subplot(3,1,2)
plot(q(pdata),ratio(offset,pdata),'k*-')
hold on
plot(q(post_data),ratio(offset,post_data),'r*-')
ylim([0,1])
xlabel('QWP angle')
ylabel('min/max ratio')
title('Power ratio')

subplot(3,1,3)
plot(q(pdata),pre_power,'k*-')
hold on
plot(q(post_data),post_power,'r*-')
legend('Pre','post')
title(sprintf('Total power, [mean,std] = [%.2f,%.2f]',mean(total_power),std(total_power)))
xlabel('QWP angle')
ylabel('Pmax+Pmin')

% subplot(4,1,4)
% title('Min/Max angle offset')
% xlabel('QWP angle')
% ylabel('Polarizer angle offset')
% suptitle('Combined')
% 
fwtext('Done!')





function corrplot(pdata,offset,ftitle,fid)
    q = @(p) p(:,1);
    
    d_th_max = wrapToPi(deg2rad(pdata(:,1) - pdata(:,3)));
    d_th_min = wrapToPi(deg2rad(pdata(:,1) - pdata(:,5)));
    total_power = pow_sum(pdata);

    sfigure(fid);
    clf;
    subplot(2,1,1)
    plot(total_power,contrast(offset,pdata),'kx')
    xlim([0,1.1*max(total_power)])
    ylim([0,1])
    xlabel('Pmin + Pmax')
    ylabel('Contrast')
    title('Contrast vs total power')
    
    subplot(2,1,2)
    plot(total_power,ratio(offset,pdata),'kx')
    xlim([0,1.1*max(total_power)])
    ylim([0,1])
    xlabel('Pmin + Pmax')
    ylabel('Ratio')
    title('Power ratio vs total power')
    
    suptitle(ftitle)
end

function bigplot(pdata,offset,ftitle,fid)

    q = @(p) p(:,1);
    
    
    d_th_max = wrapToPi(deg2rad(pdata(:,1) - pdata(:,3)));
    d_th_min = wrapToPi(deg2rad(pdata(:,1) - pdata(:,5)));
    total_power = pow_sum(pdata);
    
    sfigure(fid);
%     clf;
    subplot(4,1,1)
    plot(q(pdata),contrast(offset,pdata),'k*')
    hold on
    ylim([0,1])
    xlabel('QWP angle')
    ylabel('Contrast')
    title('Contrast')
    subplot(4,1,2)
    plot(q(pdata),ratio(offset,pdata),'k*')
    hold on
    ylim([0,1])
    xlabel('QWP angle')
    ylabel('min/max ratio')
    title('Power ratio')

    subplot(4,1,3)
    plot(q(pdata),total_power,'k*')
    hold on
    ylim([0,max(total_power)])
    title(sprintf('Total power, [mean,std] = [%.2f,%.2f]',mean(total_power),std(total_power)))
    xlabel('QWP angle')
    ylabel('Pmax+Pmin')

    subplot(4,1,4)
%     plot(q(pdata),(d_th_min),'b.')
%     hold on
%     plot(q(pdata),(d_th_max),'r.')
    plot(q(pdata),mod((d_th_max-d_th_min),pi/2)-pi/4,'k*')
    hold on
    title('Min/Max angle offset')
    xlabel('QWP angle')
    ylabel('Polarizer angle offset')
    suptitle(ftitle)
end

function P = pow_sum(pol_data)
    P = pol_data(:,3) + pol_data(:,5);
end

function C = contrast(offset,pol_data)
    C = (pol_data(:,3) - pol_data(:,5))./(pol_data(:,3) + pol_data(:,5)-2*offset);
end
function R = ratio(offset,pol_data)
    R = (pol_data(:,5)-offset)./(pol_data(:,3)-offset);
end

