
%fit probe beam
% user param

%im_raw=imread('./data/beam_profile/20180724T211918.841.png');
im_raw=imread('./data/beam_profile/20180724T211949.062.png');
pix_size=3e-6;
close('all')

%% 
beam_image=double(im_raw)/255;
beam_image=sum(beam_image,3);

thesh=min(beam_image(:))+0.2*range(beam_image(:))
area_thresh=10;
frame_fit_size_min_max=[5,1000];

frame_thresh=beam_image>thesh;
regions=regionprops('table',frame_thresh);
regions=sortrows(regions,1,'descend');
areas=regions(:,:).Area;
area_mask=areas>area_thresh;
centers = regions(:,:).Centroid;


stfig('throholding');
imagesc(beam_image)
hold on
scatter(centers(area_mask,1),centers(area_mask,2),100,'r','+')
hold off
center=regions.Centroid(1,:);
radius=regions.BoundingBox(1,3:4);
 

%%

radius=round(radius)*3;
sub_rad=max(radius); %make the bounding box square

sub_rad=max(min(sub_rad,frame_fit_size_min_max(2)),frame_fit_size_min_max(1));
sub_rad=sub_rad*[1 1];
minpos=center(1:2)-sub_rad;
minpos=max(1,minpos);
minpos=round(minpos);
maxpos=center(1:2)+sub_rad;
maxpos=min(maxpos,fliplr(size(beam_image,[1,2])));
maxpos=round(maxpos);

subframe=beam_image(minpos(2):maxpos(2),minpos(1):maxpos(1));
sub_size=size(subframe);
sub_size_m=pix_size*sub_size;
xvals_sub=linspace(-sub_size(2)/2,sub_size(2)/2,sub_size(2))*pix_size;
yvals_sub=linspace(-sub_size(1)/2,sub_size(1)/2,sub_size(1))*pix_size;

% find the mean of the image excluding this subframe
im_background=sum(beam_image(:))-sum(subframe(:));
im_background=im_background/(numel(beam_image)-numel(subframe));
subframe=subframe-im_background;
subframe=subframe./max(subframe(:));

        
%%        

x_lim=[-1,1]*30e-6;
y_lim=x_lim;
xy_factor=1e6;

font_name='cmr10';
linewidth=1.5;
font_size=12;

stfig('camera fit A');
clf
set(gca, 'FontName', font_name)
set(gca, 'FontSize', font_size)


subplot(1,2,1)
imagesc(xvals_sub*xy_factor,yvals_sub*xy_factor,subframe)
imh=gca;
pbaspect([1,1,1])
xlabel('y ($\mathrm{\mu}$m)')
ylabel('z ($\mathrm{\mu}$m)')
xlim(x_lim*xy_factor)
ylim(y_lim*xy_factor)
colormap(gca,plasma())
ch=colorbar;
ch.Label.String='Intensity (arb. u.)';
ch.FontName =font_name;
ch.Ticks=linspace(0,1,5);

subplot(1,2,2)
sh=surf(xvals_sub*xy_factor,yvals_sub*xy_factor,subframe,...
    'FaceAlpha',0.4);
ifh=gca;
shading interp
sh.EdgeColor='k';
colormap(gca,plasma())
caxis([0 0.8])
%pbaspect([1,1,1])
xlim(x_lim*xy_factor)
ylim(y_lim*xy_factor)
ifh.View= [7.8   14.04]
box on
xlabel('y ($\mathrm{\mu}$m)')
ylabel('z ($\mathrm{\mu}$m)') 
zlabel('Intensity (arb. u.)')

        
[xmesh,ymesh]=meshgrid(xvals_sub,yvals_sub);
pos_val = zeros(size(xmesh,1),size(ymesh,2),3);
pos_val(:,:,1) = xmesh;
pos_val(:,:,2) = ymesh;
pos_val(:,:,3) = subframe; 
pos_val_rows=reshape(pos_val,[],3,1);
% stfig('conversion check')
%scatter3(pos_val_rows(:,1),pos_val_rows(:,2),pos_val_rows(:,3))
%

opt = statset('TolFun',1e-10,'TolX',1e-10,...
    'MaxIter',1e4,... %1e4
    'UseParallel',1);


cof_names={'sigma x','sigma y', 'cen x','cen y', 'theta', 'offset', 'amp'};

modelfun = @(b,x) gauss_2d_rot_function(x,[b(1),b(2)],[b(3),b(4)],b(5),b(6),b(7),'waist');

fit_params_guess = [std(pos_val_rows(:,1),max(pos_val_rows(:,3),0))/5,...
                    std(pos_val_rows(:,2),max(pos_val_rows(:,3),0))/5,...
                    0,0,0,0,max(subframe(:))]; %Inital guess parameters
fitobject=fitnlm(pos_val_rows(:,1:2),pos_val_rows(:,3),modelfun,fit_params_guess,...
    'options',opt,...
    'CoefficientNames',cof_names)
fitparam=fitobject.Coefficients;
hold on
%
%subplot(1,4,3)
fit_pos_val_rows=pos_val_rows;
%fit_pos_val_rows(:,3)=modelfun(fit_params_guess,pos_val_rows(:,1:2));
fit_pos_val_rows(:,3)=predict(fitobject,pos_val_rows(:,1:2));
fit_pos_vals=reshape(fit_pos_val_rows,size(pos_val));
%imagesc(fit_pos_vals(:,:,3))
sh=surf(fit_pos_vals(:,:,1)*xy_factor,fit_pos_vals(:,:,2)*xy_factor,fit_pos_vals(:,:,3),...
    'FaceAlpha',0.8,'FaceColor',[100,200,53]./255,'edgecolor','none');
%shading interp
%colormap(viridis)



set(gcf,'Units','Pixels')
%set(gcf,'Position',[1068,355,676,453])
fig_width_px=1000;
fig_aspect_ratio=0.3; %0.67;
set(gcf,'Position',[100,355,fig_width_px,fig_width_px*fig_aspect_ratio])       
hold off

%%

fig_name='beam_profile_img';
%fig_name='pal_mean_v_acc_dyn_static';
fig_dir='./figs/thesis_figs';
export_fig(fullfile(fig_dir,strcat(fig_name,'.svg')))
export_fig(fullfile(fig_dir,strcat(fig_name,'.pdf')))

%%


%colors_main=[[98,136,63];[95,109,187];[180,72,117]]./255;
colors_main=[[40,136,40];[95,109,187];[218, 129, 57]]./255;

hsv=colorspace('RGB->HSV',colors_main(:,:));
hsv(:,2)=hsv(:,2);
colors_shaded=colorspace('HSV->RGB',hsv);

stfig('camera fit B')
clf
set(gca, 'FontName', font_name)
set(gca, 'FontSize', font_size)


% [fit_params_out,resnorm,residual,exitflag]=lsqcurvefit(,fit_params_guess,...
%     pos_val,subframe,options);

    
%subplot(2,4,6)
%imagesc(xvals_sub,yvals_sub,D2GaussFunctionRot(fit_params_out,pos_data)) %plot fit


subplot(1,2,1)
%x line profile
x1d_plot_vals=linspace(min(x_lim),max(x_lim),1e3);
xysamples=[];
xysamples(:,1)=x1d_plot_vals;
xysamples(:,2)=x1d_plot_vals*0+fitobject.Coefficients{'cen y','Estimate'}+-2e-6;
[~,ci]=predict(fitobject,xysamples,'Alpha',1-erf(1/sqrt(2)),'Prediction','curve');
[prediction,oi]=predict(fitobject,xysamples,'Alpha',1-erf(1/sqrt(2)),'Prediction','observation');
measured_vals=interp2(pos_val(:,:,1),pos_val(:,:,2),pos_val(:,:,3),xysamples(:,1),xysamples(:,2));


ci_up=oi(:,1);
ci_down=oi(:,2);
patch([x1d_plot_vals, fliplr(x1d_plot_vals)]*xy_factor, ...
    [ci_up', fliplr(ci_down')], colors_shaded(1,:),'EdgeColor','none',...
    'FaceAlpha',0.3)  %[1,1,1]*0.80
hold on
% ci_up=ci(:,1);
% ci_down=ci(:,2);
% patch([x1d_plot_vals, fliplr(x1d_plot_vals)]*xy_factor, ...
%     [ci_up', fliplr(ci_down')], colors_shaded(3,:),'EdgeColor','none',...
%     'FaceAlpha',0.5)  %[1,1,1]*0.80

plot(x1d_plot_vals*xy_factor,measured_vals,'r','LineWidth',linewidth,'Color',colors_main(3,:))
plot(x1d_plot_vals*xy_factor,prediction,'--','LineWidth',linewidth,'Color',colors_main(1,:))
hold off
xlim(x_lim*xy_factor)
xlabel('y ($\mathrm{\mu}$m)')
ylabel('Intensity (arb. u.)')
box on
legend('','camera','fit')
        
subplot(1,2,2)
y1d_plot_vals=linspace(min(y_lim),max(y_lim),1e3);
xysamples=[];
xysamples(:,1)=x1d_plot_vals*0+fitobject.Coefficients{'cen x','Estimate'}+0e-6;
xysamples(:,2)=y1d_plot_vals;

[~,ci]=predict(fitobject,xysamples,'Alpha',1-erf(1/sqrt(2)),'Prediction','curve');
[prediction,oi]=predict(fitobject,xysamples,'Alpha',1-erf(1/sqrt(2)),'Prediction','observation');
measured_vals=interp2(pos_val(:,:,1),pos_val(:,:,2),pos_val(:,:,3),xysamples(:,1),xysamples(:,2));
ci_up=oi(:,1);
ci_down=oi(:,2);
patch([x1d_plot_vals, fliplr(x1d_plot_vals)]*xy_factor, ...
    [ci_up', fliplr(ci_down')], colors_shaded(1,:),'EdgeColor','none',...
    'FaceAlpha',0.3)  %[1,1,1]*0.80
hold on
plot(x1d_plot_vals*xy_factor,measured_vals,'r','LineWidth',linewidth,'Color',colors_main(3,:))
plot(x1d_plot_vals*xy_factor,prediction,'--','LineWidth',linewidth,'Color',colors_main(1,:))
hold off
xlim(y_lim*xy_factor)
xlabel('z ($\mathrm{\mu}$m)')
ylabel('Intensity (arb. u.)')
box on
legend('','camera','fit')


 
radi=[fitobject.Coefficients{'sigma x','Estimate'},fitobject.Coefficients{'sigma y','Estimate'}];
radi=sort(radi);
fprintf('beam radi %s um , rms value %.3f um \n',sprintf('%.3f,',radi*1e6),rms(radi)*1e6)

set(gcf,'Units','Pixels')
%set(gcf,'Position',[1068,355,676,453])
fig_width_px=1000;
fig_aspect_ratio=0.3; %0.67;
set(gcf,'Position',[100,355,fig_width_px,fig_width_px*fig_aspect_ratio])       
hold off

%%
%grid on
fig_name='beam_profile_line';
%fig_name='pal_mean_v_acc_dyn_static';
fig_dir='./figs/thesis_figs';
export_fig(fullfile(fig_dir,strcat(fig_name,'.svg')))
export_fig(fullfile(fig_dir,strcat(fig_name,'.pdf')))
