%beam profiler

%user settings
pix_size=3e-6;

%end user settings
%%
objects = imaqfind;
if ~isequal(objects,[])
    delete(objects);
end

info=imaqhwinfo('winvideo');
formats=info.DeviceInfo.SupportedFormats;
    
vid = videoinput('winvideo', 1,'MJPG_1920x1080');
vid.FramesPerTrigger = 1;
triggerconfig(vid,'manual')
addpath('Colormaps')

%%
preview(vid);
%%
src = getselectedsource(vid);
vid.FramesPerTrigger = 1;

src.BacklightCompensation='off';
src.ExposureMode='manual';
src.Exposure=-13; %[-13 -1] %t=2^exposure
src.Contrast=0; %[0 64]
src.Sharpness=0;
src.Gamma=100; %vout+A vin ^ 100/gamma [72 500]
src.Gain=50; %vout+A vin ^ 100/gamma [72 500]
src.Brightness=-64; %[-64 64]
src.Hue=0;
src.WhiteBalanceMode='manual';
pause(0.1)
%%

exp_prop=propinfo(src,'Exposure');
bright_prop=propinfo(src,'Brightness');

figure(2)
bright_val=20;
bright_lim=bright_prop.ConstraintValue;
exp_lim=exp_prop.ConstraintValue;
target_maxval=0.5;
src.Exposure=-4; 
blur=5;
thesh=5*max(max(fspecial('gaussian',blur*5,blur)));
area_thresh=10;

pidstate.integrator=

while true
    
    src.gain=bright_val;
    frame_sd = double(getsnapshot(vid))./double(intmax('uint8')); %scaled and in double format
    frame_sd=double(squeeze(frame_sd(:,:,1)));
    frame_sd=flipud(frame_sd');
    frame_sd=frame_sd-imgaussfilt(frame_sd,300,'Padding','circular');
    stat=[max(frame_sd(:)),mean(frame_sd(:)),std(single(frame_sd(:)))];
    frame_smooth=imgaussfilt(frame_sd,blur,'Padding','circular');

    frame_thresh=frame_smooth>thesh;
    regions=regionprops('table',frame_thresh);
    regions=sortrows(regions,1,'descend');
    areas=regions(:,:).Area;
    area_mask=areas>area_thresh;
    centers = regions(:,:).Centroid;
    
    
    subplot(2,4,1)
    imagesc(frame_sd)
    colormap(gca,inferno())
    caxis([0 1])
    pbaspect([1,1,1])

    subplot(2,4,2)
    imagesc(frame_smooth)
    colormap(gca,inferno())
    caxis([0 1])
    pbaspect([1,1,1])

    subplot(2,4,3)
    imagesc(frame_thresh)
    colormap(inferno())  
    fit_params_out=NaN(7,1); %initalize
    if ~isequal(centers,[])
        hold on
        scatter(centers(area_mask,1),centers(area_mask,2),100,'r','+')
        hold off
        center=round(regions.Centroid(1,:));
        radius=round(regions.BoundingBox(1,3:4));
        sub_rad=max(radius*2); %make the bounding box square
        sub_rad=max(min(sub_rad,100),10);
        sub_rad=sub_rad*[1 1];
        minpos=center-sub_rad;
        minpos=max(1,minpos);
        maxpos=center+sub_rad;
        maxpos=min(maxpos,fliplr(size(frame_sd)));
        
        
        subframe=frame_sd(minpos(2):maxpos(2),minpos(1):maxpos(1));
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
        fit_params_guess = [1,0,radius(1)/2,0,radius(2)/2,0,0]; %Inital guess parameters
        fit_params_lb = [0,-sub_size(1)/2,0,-sub_size(2)/2,0,-pi/4,-0.5]; 
        fit_params_ub = [5,sub_size(1)/2,(sub_size(1)/2)^2,sub_size(1)/2,(sub_size(2)/2)^2,pi/4,0.5];
        options=optimoptions('lsqcurvefit','FunctionTolerance',1e-4,'Display','off');
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
        
    else
        subplot(2,4,4)
        imagesc(zeros(10))
        caxis([0 1])
        
        subplot(2,4,5)
        surf(zeros(10))
    end
    
    fprintf('max amp %01.3f, b %+02i e %+02i fit wx %01.3f wy %01.3f µm \n',...
        stat(1),src.Brightness,src.Exposure,...
        fit_params_out(3)*pix_size*1e6,fit_params_out(5)*pix_size*1e6)
    
    [ctr_output,pidstate]=pid(pidstate);
    
    bright_val=int32(round(double(bright_val)+(target_maxval-stat(1))*50));
    bright_val=min(max(bright_lim(1),bright_val),bright_lim(2));
    
    %if stat(1)>target_maxval && bright>blim(1)
    %    bright=bright-1;
    %elseif stat(1)<target_maxval && bright<blim(2)
    %    bright=bright+1;
    %end
   
    pause(1e-6)
end
fprintf('\n')







%%
 %test a auto exposure system
figure(1)

gain_list=linspace(min(gain_lim),max(gain_lim),50);
frame_struct=[];
cam.Brightness=0;
cam.Exposure=-13;
cam.Brightness=-40;
gain=10;
pidstate.integrator=(gain-min(gain_lim))/range(gain_lim);
pidstate.setpt=0.85;
pidstate.k_int=-1e-6;
pidstate.k_prop=1e-1;
pidstate.outlims=[5/100 1];
pidstate.aw_thresh_range=0.05; %how far away from the edge aw starts (full range 0-1)
pidstate.int_lim=3;
pidstate.slew_lim=1e-5;
ii=1;
maxval=0;
exposure=cam.Exposure;

for n=1:100
    cam.Gain=gain;
    fprintf('exp %2i, gain %2.1f, maxval %2.3f int %2.3f \n',[cam.Exposure,cam.Gain,maxval,pidstate.integrator])
    frame_struct.setings(ii,:)=[cam.Exposure cam.Gain cam.Brightness];
    pause(0.001)
    frame = double(snapshot(cam))/(2^8-1);
    frame=squeeze(frame(:,:,1));
    [N,edges] = histcounts(frame(:),linspace(0,1,1e5));
    centers=edges(1:end-1)+diff(edges)/2;
    frame_struct.hist(ii,:)=[N,centers];
    maxval=max(frame(:));
    frame_struct.stat(ii,:)=[min(frame(:)),max(frame(:)),mean(frame(:)),std(frame(:))];
    pidstate.meas=maxval;
    pidstate=pid_loop(pidstate);
    gain=pidstate.ctr_output*range(gain_lim)+min(gain_lim);
    ii=ii+1;
    if pidstate.aw<0.4 && exposure~=exposure_lim((pidstate.integrator>0.5)*1+1)
        exposure=cam.Exposure+sign(pidstate.integrator);
        cam.Exposure=exposure;
        gain=50;
        pidstate.integrator=(gain-min(gain_lim))/range(gain_lim);
    end
end









%%
delete(vid);