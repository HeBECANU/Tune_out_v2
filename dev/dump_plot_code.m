file_num=1
fignum=1
anal_opts.appr_freq_guess=[52,47.6,40];


velocity=const.g0*anal_opts.fall_time
trange_single_sec=range(import_opts.txylim(1,:));
trange_single_mm=trange_single_sec*velocity;
xrange_single_mm=range(histplot.xlim);
histplot.binst=round(histplot.binsx*trange_single_mm/xrange_single_mm)+1;
XEdges=linspace(histplot.xlim(1),histplot.xlim(2),histplot.binsx);
TEdges=linspace(import_opts.txylim(1,1),import_opts.txylim(1,2),histplot.binst);
TEdges=linspace(import_opts.txylim(1,1),import_opts.txylim(1,2),histplot.binsx);
%counts_all=cat(1,data.txy{:});
counts_all= data.mcp_tdc.counts_txy{file_num};

idx=[2,3,2];
[counts,centers]=hist3(counts_all(:,[1,idx(histplot.dimesion)]),'edges',{TEdges,XEdges}); %[min(TEdges),min(YEdges)]
%fprintf('maxcounts in a bin %f \n',max(max(counts)))

%fprintf('max density %f \n',max_density_presmooth)
    %counts=counts/max(counts(:));
if  ~histplot.blur==0
    filt_ratio=(trange_single_mm/xrange_single_mm)*(histplot.binsx/histplot.binst);
    filt_size=round(100*histplot.blur);
    filt_size=[round(filt_ratio*filt_size) filt_size];
    filter=Gaussian_filter2d(filt_size,[filt_ratio*histplot.blur histplot.blur]);
    counts_smooth=imfilter(counts, filter, 'same');
else
    counts_smooth=counts;
end

fig=figure(fignum);
subplot(3,1,1)
%imagesc seems to plot the wrong way round so we transpose here
imagesc(centers{1},1e3*centers{2},counts_smooth');
yl=ylim;
xl=xlim;
%can add line on top
hold on
plot(squeeze(anal_out.mean_pos_window(file_num,:,3)),1e3*squeeze(anal_out.mean_pos_window(file_num,:,2+idx(histplot.dimesion))),'rx','MarkerSize',12)
plot(squeeze(anal_out.mean_pos_window(file_num,:,2)),1e3*squeeze(anal_out.mean_pos_window(file_num,:,2+idx(histplot.dimesion))),'wx','MarkerSize',12)
hold off
title(sprintf('file %03i',file_num))
ylabel('Y Pos (mm)')
legend('center','win cen')
xlabel('Time (s)')
set(gca,'Ydir','normal')
set(gcf,'Color',[1 1 1]);
%pbaspect([trange_single_mm xrange_single_mm 1])
%xticklabels({'','','','',''})
cm_treshold=5;
cm_shift=10;
cm_var=inferno();%viridis();
cm_bg=[repmat(cm_var(1,:),cm_treshold,1) ; cm_var(cm_shift:end,:)];
colormap(cm_var)
%yyaxis right
%ylim(yl*fall_time)
%ylabel('X Position(mm)')
%set(gcf, 'Units', 'Pixels', 'OuterPosition', [100, 100, 500, 300])

subplot(3,1,2)
xvalues=1e3*squeeze(anal_out.mean_pos_window(file_num,:,4))';
xvalues=xvalues-mean(xvalues);
yvalues=1e3*squeeze(anal_out.mean_pos_window(file_num,:,5))';
yvalues=yvalues-mean(yvalues);
zvalues=squeeze(anal_out.mean_pos_window(file_num,:,3))-squeeze(anal_out.mean_pos_window(file_num,:,2));
zvalues=zvalues'*velocity*1e3;
zvalues=zvalues-mean(zvalues);

sqrtn=sqrt(anal_out.mean_pos_window(file_num,:,9)); %find the statistical uncert in a single shot
sqrtn=0.1*sqrtn;
tvalues=(1:numel(anal_out.mean_pos_window(file_num,:,1)))'*anal_opts.pulsedt;
xerr=1e3*squeeze(anal_out.mean_pos_window(file_num,:,7))'./sqrtn;
yerr=1e3*squeeze(anal_out.mean_pos_window(file_num,:,8))'./sqrtn;
zerr=1e3*squeeze(anal_out.mean_pos_window(file_num,:,6))'*velocity./sqrtn;
xyzerr=[xerr,yerr,zerr];

plot(tvalues,xvalues,'kx-')
hold on
plot(tvalues,yvalues,'rx-')
plot(tvalues,zvalues,'bx-')
%plot(squeeze(mean_pos_window(1,1,:,1)),...
%    squeeze(mean_pos_window(1,1,:,3)-mean_pos_window(1,1,:,1)+mean_pos_window(1,1,1,1))*velocity*1e3,'kx')
hold off
ylabel('X Pos (mm)')
xlabel('Time (s)')
set(gca,'Ydir','normal')
set(gcf,'Color',[1 1 1]);
legend('x','y','z')
pause(0.1)



subplot(3,1,3)
%modelfun = @(b,x) b(1)*sin(b(2)*x(:,1)*pi+b(3)*pi)+b(4)+b(5)*x(:,1)+x(:,2)*b(6)+x(:,3)*b(7);
%beta0=[15, 30, -0.5, -2, 0, 0, 0];
%cof_names={'amp','freq','phase','offset','grad','ycpl','zcpl'}

modelfun = @(b,x) exp(-x(:,1)./max(0,b(7))).*b(1).*sin(b(2)*x(:,1)*pi*2+b(3)*pi*2)+b(4)+b(8)*x(:,1)+b(5)*x(:,2)*0+b(6)*x(:,3)*0;
beta0=[5, anal_opts.appr_freq_guess(histplot.dimesion), 1, 1,0,0,0.5,0.1];
cof_names={'amp','freq','phase','offset','ycpl','zcpl','damp','grad'};

opt = statset('TolFun',1e-10,...
    'TolX',1e-10,...
    'MaxIter',1e4,...
    'UseParallel',1);
txyz=[tvalues,xvalues,yvalues,zvalues];
idx=1:4;
idx(histplot.dimesion+1)=[];
predictor=txyz(:,idx);
%predictor=[tvalues,xvalues,zvalues];
fitobject=fitnlm(predictor,txyz(:,histplot.dimesion+1),modelfun,beta0,'options',opt,...
    'CoefficientNames',cof_names)

fitparam=fitobject.Coefficients;
anal_out.fitparams{file_num}=fitparam;
%limiting frequnecy prediction from http://adsabs.harvard.edu/full/1999DSSN...13...28M
meanwidth=sqrt(mean(squeeze(anal_out.mean_pos_window(file_num,:,7)).^2))*1e3;
%meanwidth=mean(xyz_unc(:,dimension))

frequnclim=sqrt(6/sum(anal_out.mean_pos_window(file_num,:,9)))*...
    (1/(pi*range(tvalues)))*...
    (meanwidth/fitparam{2,1});
fprintf('sampling limit %2.3g Hz, fint unc %2.3g Hz, ratio %2.3g \n',[frequnclim,fitparam{2,2},fitparam{2,2}/frequnclim])

anal_out.fit_sample_limit{file_num}=[frequnclim,fitparam{2,2},fitparam{2,2}/frequnclim];


%could include a cross freq term for the z and y components

deltphase=0.05*2*pi;
fprintf('xy angle %2.5f° \n',asin(fitparam{5,1})*(180/pi))
fprintf('xz angle %2.5f° \n',asin(fitparam{6,1})*(180/pi))
fprintf('time untill phase unc reaches %2.3f°  %2.3fcyc= %2.3f \n',...
    deltphase*(180/pi),deltphase/(2*pi),(deltphase-fitparam{3,2})/fitparam{2,2})

tplotvalues=linspace(min(tvalues),max(tvalues),1e5)';
predictorplot=[tplotvalues,...
               interp1(predictor(:,1),predictor(:,2),tplotvalues),...
               interp1(predictor(:,1),predictor(:,3),tplotvalues)];
[prediction,ci]=predict(fitobject,predictorplot);
plot(predictorplot(:,1),prediction,'-','LineWidth',1.5,'Color',[0.5 0.5 0.5])
ax = gca;
set(ax, {'XColor', 'YColor'}, {'k', 'k'});
hold on
plot(predictorplot(:,1),ci(:,1),'-','LineWidth',1.5,'Color','k')
plot(predictorplot(:,1),ci(:,2),'-','LineWidth',1.5,'Color','k')
errorbar(predictor(:,1),txyz(:,histplot.dimesion+1),xyzerr(:,histplot.dimesion),'k.','MarkerSize',10,'CapSize',0,'LineWidth',1,'Color','r') 
set(gcf,'Color',[1 1 1]);
ylabel('X(mm)')
xlabel('Time (s)')
hold off
ax = gca;
set(ax, {'XColor', 'YColor'}, {'k', 'k'});
set(gca,'linewidth',1.0)
%set(gcf, 'Units', 'Pixels', 'OuterPosition', [100, 100, 500, 300])
%set(gca,'ylim',[-20,20])
%yl=ylim;
%yticks(linspace(yl(1),yl(2),6))
%yticks([-8 -4 0 4 8])
pause(0.01)
