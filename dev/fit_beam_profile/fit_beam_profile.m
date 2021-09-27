
%fit probe beam
beam_image=imread('./data/20180724T211918.841.png');
beam_image=double(beam_image)/255;
beam_image=sum(beam_image,3);

thesh=min(beam_image(:))+0.2*range(beam_image(:))
area_thresh=10;
frame_fit_size_min_max=[10,1000];
pix_size=3e-6;

frame_thresh=beam_image>thesh;
regions=regionprops('table',frame_thresh);
regions=sortrows(regions,1,'descend');
areas=regions(:,:).Area;
area_mask=areas>area_thresh;
centers = regions(:,:).Centroid;


stfig('throholding')
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
        
%%        
subframe=beam_image(minpos(2):maxpos(2),minpos(1):maxpos(1));
sub_size=size(subframe);
sub_size_m=pix_size*sub_size;

xvals_sub=linspace(-sub_size(2)/2,sub_size(2)/2,sub_size(2));
yvals_sub=linspace(-sub_size(1)/2,sub_size(1)/2,sub_size(1));

subplot(2,4,4)
imagesc(xvals_sub,yvals_sub,subframe)
pbaspect([1,1,1])

subplot(2,4,5)
surf(subframe)
colormap(gca,viridis())
caxis([0 1])
shading interp
%pbaspect([1,1,1])
        
        
[xmesh,ymesh]=meshgrid(xvals_sub,yvals_sub);
pos_data = zeros(size(xmesh,1),size(ymesh,2),2);
pos_data(:,:,1) = xmesh;
pos_data(:,:,2) = ymesh;

%%
options=optimoptions('lsqcurvefit','FunctionTolerance',1e-2,'Display','off','MaxFunctionEvaluations',10000);
fit_params_guess = [1,0,radius(1)/2,0,radius(2)/2,0,0]; %Inital guess parameters
fit_params_lb = [0,-sub_size(1)/2,0,-sub_size(2)/2,0,-pi/4,-0.5]; 
fit_params_ub = [5,sub_size(1)/2,(sub_size(1)/2)^2,sub_size(1)/2,(sub_size(2)/2)^2,pi/4,0.5];
[fit_params_out,resnorm,residual,exitflag]=lsqcurvefit(@D2GaussFunctionRot,fit_params_guess,...
    pos_data,subframe,fit_params_lb,fit_params_ub,options);

    
subplot(2,4,6)
imagesc(xvals_sub,yvals_sub,D2GaussFunctionRot(fit_params_out,pos_data)) %plot fit


subplot(2,4,7)
%x line profile
x1d_plot_vals=linspace(min(xvals_sub),max(xvals_sub),1e3);
xysamples=[];
xysamples(:,:,1)=x1d_plot_vals;
xysamples(:,:,2)=x1d_plot_vals*0;
fitvals=D2GaussFunctionRot(fit_params_out,xysamples);
measured_vals=interp2(xvals_sub,yvals_sub,subframe,xysamples(:,:,1),xysamples(:,:,2));
plot(x1d_plot_vals,fitvals)
hold on
plot(x1d_plot_vals,measured_vals,'r')
hold off
        
subplot(2,4,8)
y1d_plot_vals=linspace(min(yvals_sub),max(yvals_sub),1e3);
xysamples=[];
xysamples(:,:,1)=y1d_plot_vals*0;
xysamples(:,:,2)=y1d_plot_vals;
fitvals=D2GaussFunctionRot(fit_params_out,xysamples);
measured_vals=interp2(xvals_sub,yvals_sub,subframe,xysamples(:,:,1),xysamples(:,:,2));
plot(y1d_plot_vals,fitvals)
hold on
plot(y1d_plot_vals,measured_vals,'r')
hold off

cen_fit=fit_params_out([2,4])+center(1:2); %absolute position on ccd in px
 
radi=[fit_params_out(3),fit_params_out(5)];
radi=sort(radi)*pix_size;
fprintf('beam radi %s um , rms value %.3f um \n',sprintf('%.3f,',radi*1e6),rms(radi)*1e6)

            