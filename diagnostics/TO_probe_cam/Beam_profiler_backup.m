%beam profiler

%user settings
pix_size=3e-6;

%end user settings
%%
addpath('Colormaps')
clear('cam')
cam = webcam('HD USB Camera');
resolutions=cam.AvailableResolutions;
cam.Resolution=resolutions{1};
cam.BacklightCompensation=0;
cam.ExposureMode='manual';
cam.Exposure=-1; %[-13 -1] %t=2^exposure
cam.Contrast=0; %[0 64]
cam.Sharpness=0;
cam.Gamma=100; %vout+A vin ^ 100/gamma [72 500]
cam.Brightness=10; %[-64 64]
cam.Hue=0;
cam.WhiteBalanceMode='manual';
%preview(cam)


%webcam has no method for allowed values so i these will be hard coded for now
brightness_lim=[-64 64];
exposure_lim=[-13 -2];
gamma_lim=[0 100];
gain_lim=[0 100];

%%
%preview(cam);
%%


%%
pause(0.1)
blur=5;
fit_every=4; %fit every n itterations
thesh=0.4;%5*max(max(fspecial('gaussian',blur*5,blur)));
area_thresh=10;
cam.Exposure=-13;
cam.Brightness=-40;
cam.Sharpness=0;
gain=12; %inital gain

pidstate.integrator=(gain-min(gain_lim))/range(gain_lim);
pidstate.setpt=0.50; %max val set pt
pidstate.k_int=-2e-7;
pidstate.k_prop=1e-5;
pidstate.outlims=[12/100 1];
pidstate.aw_thresh_range=0.05; %how far away from the edge aw starts (full range 0-1)
pidstate.int_lim=3;
pidstate.slew_lim=1e-5;
center_mem=zeros(20,2); %needs to be even length
radius_mem=center_mem+10;
center_mem=center_mem+500;
ii=1;
maxval=0;
logist=@(x) 1./(1+exp(-x));


fit_params_out=NaN(7,1);
fit_params_guess = [1,0,radius(1)/2,0,radius(2)/2,0,0]; %Inital guess parameters
fit_params_lb = [0,-sub_size(1)/2,0,-sub_size(2)/2,0,-pi/4,-0.5]; 
fit_params_ub = [5,sub_size(1)/2,(sub_size(1)/2)^2,sub_size(1)/2,(sub_size(2)/2)^2,pi/4,0.5];
options=optimoptions('lsqcurvefit','FunctionTolerance',1e-4,'Display','off');
tic;
figure(2)

while true
    frame_sd = double(snapshot(cam))/(2^8-1); %scaled and in double format
    frame_sd=squeeze(frame_sd(:,:,1));
    frame_sd=flipud(frame_sd');
    maxval=max(frame_sd(:));
    stat=[maxval,mean(frame_sd(:)),std(single(frame_sd(:)))];
    %frame_sd=frame_sd-imgaussfilt(frame_sd,200,'Padding','circular');
    frame_smooth=frame_sd;
    %frame_smooth=imgaussfilt(frame_sd,blur,'Padding','circular');
    pidstate.meas=maxval;
    pidstate=pid_loop(pidstate);
    gain=pidstate.ctr_output*range(gain_lim)+min(gain_lim);
    if pidstate.aw<0.1 && exposure~=exposure_lim((pidstate.integrator>0.5)*1+1)
        exposure=cam.Exposure+sign(pidstate.integrator-0.5);
        gain=10;
        pidstate.integrator=(gain-min(gain_lim))/range(gain_lim);
        pause(0.01)
        cam.Exposure=exposure;
    end
    cam.Gain=gain;
    
    
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
     %initalize
    if ~isequal(centers,[]) 
        hold on
        scatter(centers(area_mask,1),centers(area_mask,2),100,'r','+')
        hold off
        center=regions.Centroid(1,:);
        radius=regions.BoundingBox(1,3:4);
        
%         center_mem=circshift(center_mem,1,1);
         radius_mem=circshift(radius_mem,1,1);
%         center_mem(1,:)=center;
         radius_mem(1,:)=radius;
%         fprintf('%f\n',std(center_mem(:,1)))
%         %nonlinear filter, change the smoothing time based on the std of the cen hist
%         %larger the std the shorter the filtering time
%         %could add in a hysterisis term where values only update if they change more than a ceritan
%         %thershold
%         fac=0.1;
%         winsize=@(x) 10*logist((x-4)*2);
%         win=[gausswin(size(center_mem,1)*2,winsize(nanstd(center_mem(1:10,1)))),...
%             gausswin(size(center_mem,1)*2,winsize(nanstd(center_mem(1:10,2))))];
%         win=win(1+size(win,1)/2:end,:);
%         win=[win(:,1)./nansum(win(:,1)),win(:,2)./nansum(win(:,2))];
%         cen_filt=nansum(center_mem.*win);
        %center=round(cen_filt(1,:));
        radius=round(mean(radius_mem));
        center=round(center);
        radius=round(radius);
        sub_rad=max(radius*3); %make the bounding box square
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

        if mod(ii,fit_every)==0
            [fit_params_out,resnorm,residual,exitflag]=lsqcurvefit(@D2GaussFunctionRot,fit_params_guess,...
                pos_data,subframe,fit_params_lb,fit_params_ub,options);
        end
        
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
        
        cen_fit=fit_params_out([2,4])+center; %absolute position on ccd in px
    else
        fit_params_out=NaN(7,1);
        subplot(2,4,4)
        imagesc(zeros(10))
        caxis([0 1])
        
        subplot(2,4,5)
        surf(zeros(10))
        cen_fit=[NaN NaN];
    end
    pause(1e-6)
    ii=ii+1;
    radi=[fit_params_out(3),fit_params_out(5)];
    radi=sort(radi)*pix_size*1e6;
    
    fprintf('max amp %04.3f, g %05.2f e -%02i fit wx %05.2f wy %05.2f µm x %05.2f y %05.2f px fr %03.1f \n',...
        maxval,cam.Gain,-cam.Exposure,...
        radi(1),radi(2),cen_fit(1),cen_fit(2),1/toc)
    tic
end
fprintf('\n')






%%
delete(vid);
