%beam profiler


%end user settings
%%
clc; clear;
fclose('all')
objects = imaqfind;
delete(objects);


%%
this_folder = fileparts(which(mfilename));
% Add that folder plus all subfolders to the path.
addpath(genpath(this_folder));

info=imaqhwinfo('winvideo');
formats=info.DeviceInfo.SupportedFormats;
vid = videoinput('winvideo', 1,'MJPG_1920x1080');
triggerconfig(vid,'immediate')
vid.TriggerRepeat=inf;
vid.FramesPerTrigger = 1;
preview(vid) %critical to have this for good performance, it seems matlab makes an iternal buffer when preview is used
%which does not block the setting of exposure

cam = getselectedsource(vid);
cam.BacklightCompensation='off';
cam.ExposureMode='manual';
cam.Exposure=-13; %[-13 -1] %t=2^exposure
cam.Contrast=0; %[0 64]
cam.Sharpness=0;
cam.Gamma=100; %vout+A vin ^ 100/gamma [72 500]
cam.Brightness=10; %[-64 64]
cam.Hue=0;
cam.WhiteBalanceMode='manual';




%webcam has no method for allowed values so i these will be hard coded for now
property=propinfo(cam,'Exposure');
exposure_lim=double(property.ConstraintValue);
property=propinfo(cam,'Brightness');
brightness_lim=double(property.ConstraintValue);
property=propinfo(cam,'Exposure');
exposure_lim=property.ConstraintValue;
property=propinfo(cam,'Gamma');
gamma_lim=property.ConstraintValue;
property=propinfo(cam,'Gain');
gain_lim=double(property.ConstraintValue);


%%

%START USER PARAM
pix_size=3e-6;
blur=5;
fit_every=1; %fit every n itterations
write_log_every=1; %write every n seconds
save_image_every=1;
thesh=0.3;%5*max(max(fspecial('gaussian',blur*5,blur)));
area_thresh=10;
exposure=-12;
cam.Brightness=-40;
cam.Sharpness=0;
frame_fit_size_min_max=[10,80];
gain=8; %inital gain
%END USER PARAM

pidstate=[];
%set the inital output
pidstate.ctr_output=exposure;
cam.Exposure=pidstate.ctr_output;
%set the pid params
pidstate.initalize=1;
pidstate.setpt=0.60;
pidstate.verbose=0;
pidstate.k_int=-1e2;
pidstate.k_prop=-3e-2;
pidstate.outlims=[5,95];
pidstate.aw_thresh_range=0.05; %how far away from the edge AW starts 
%pidstate.int_lim=3;
pidstate.slew_lim=20;
center_mem=zeros(20,2); %needs to be even length
radius_mem=center_mem+10;
center_mem=center_mem+500;
ii=1;
maxval=0;
logist=@(x) 1./(1+exp(-x));

time_start_posix=posixtime(datetime('now'));
time_last_image_save=time_start_posix;
time_last_log=time_start_posix;
fit_params_out=NaN(7,1);
logfile= fopen('./logs/beam_profiler_log.txt','a');
options=optimoptions('lsqcurvefit','FunctionTolerance',1e-2,'Display','off','MaxFunctionEvaluations',10);
tic;
figure(2)

while true
    
    frame_raw =getsnapshot(vid);
    time_now_posix=posixtime(datetime('now'));
    frame_sd = double(frame_raw)/(2^8-1); %scaled and in double format
    frame_sd=squeeze(frame_sd(:,:,1));
    frame_sd=flipud(frame_sd');
    maxval=max(frame_sd(:));
    stat=[maxval,mean(frame_sd(:)),std(single(frame_sd(:)))];
    %frame_sd=frame_sd-imgaussfilt(frame_sd,200,'Padding','circular');
    frame_smooth=frame_sd;
    %frame_smooth=imgaussfilt(frame_sd,blur,'Padding','circular');
    pidstate.meas=maxval;
    pidstate=pid_loop(pidstate);
    if pidstate.aw<0.4 && exposure~=exposure_lim((pidstate.integrator>50)*1+1) 
        exposure=cam.Exposure+sign(pidstate.integrator-0.5);
        if pidstate.integrator-0.5>1
            pidstate.ctr_output=20;
        else
            pidstate.ctr_output=90;
        end
        pidstate.initalize=true;
        cam.Exposure=exposure;
    end
    cam.Gain=pidstate.ctr_output;
    
    
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
        sub_rad=max(radius*1); %make the bounding box square
        
        sub_rad=max(min(sub_rad,frame_fit_size_min_max(2)),frame_fit_size_min_max(1));
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
            fit_params_guess = [1,0,radius(1)/2,0,radius(2)/2,0,0]; %Inital guess parameters
            fit_params_lb = [0,-sub_size(1)/2,0,-sub_size(2)/2,0,-pi/4,-0.5]; 
            fit_params_ub = [5,sub_size(1)/2,(sub_size(1)/2)^2,sub_size(1)/2,(sub_size(2)/2)^2,pi/4,0.5];
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
    if write_log_every==inf
        write_log=false;
    else
        write_log=time_now_posix>time_last_log+write_log_every;
    end
    if save_image_every==inf
        write_image=false;
    else
        write_image=time_now_posix>time_last_image_save+save_image_every;
    end
    %
    %%
    if write_log || write_image
        
        %com: this should be changed to json format
        %write to a csv file
        %all times in canberra local time AEST
        %row feilds: posix_time(ms),hr_time(iso),max
        %amp,gain,exposure,brightness,rmin,rmax,cenx,ceny
         nowdt=datetime('now');
         logstr=sprintf('%.3f,%s,%05.3f,%03u,%03i,%2i,%07.3f,%07.3f,%08.3f,%08.3f\n',...
             posixtime(nowdt),datestr(nowdt,'yyyy-mm-ddTHH:MM:SS.FFF'),maxval,...
             cam.Gain,cam.Exposure,cam.Brightness,radi(1),radi(2),cen_fit(1),cen_fit(2));
        fprintf(logfile,logstr);
        time_last_log=time_now_posix;
    end    
    if write_image
        
         imwrite(frame_raw,['./logs/',datestr(datetime('now'),'yyyymmddTHHMMSS.FFF'),'.png'],...
             'Mode','lossless','Compression','deflate','Comment',logstr,'BitDepth',8)
        time_last_image_save=time_now_posix;
        
    end
    %%
    fprintf('max%05.3f,g%03u,e%03i,b%2i,fitwx%06.2f,wy%06.2fµm,x%07.2f,y%07.2fpx,fr%03.1f \n',...
        maxval,cam.Gain,cam.Exposure,cam.Brightness,...
        radi(1),radi(2),cen_fit(1),cen_fit(2),1/toc)
    tic
end
fprintf('\n')

%%

tic
times=zeros(200,1);
for n=1:size(times,1)
    cam.Gain=mod(n,10)*10;
    frame = getsnapshot(vid);
    imagesc(frame)
    times(n)=toc;
    tic
end
mean(1./times)




%%
fclose(logfile)
delete(vid);
clear vid;
