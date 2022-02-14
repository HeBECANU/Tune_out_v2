%some scripts to play with http://www.elpcctv.com/full-hd-2mp-usb-camera-black-and-white-module-p-231.html
%eventualy want to build into a program to be used for beam profiling
%find the peak region, add lines to main plot
%subplot of this peak region with gaussian fit (line/2d)
%automatic exposure control and linearization to give a measure of power
%hilight saturated regions
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
preview(cam)


%webcam has no method for allowed values so i these will be hard coded for now
brightness_lim=[-64 64];
exposure_lim=[-13 -2];
gamma_lim=[0 100];
gain_lim=[0 100];


 %%
frame = snapshot(cam);
frame=squeeze(frame(:,:,1));
imagesc(flipud(frame'))
%colormap('veridis')
colormap(inferno())
caxis([0 intmax('uint8')])
pbaspect([1,1,1])
size(frame)
max(max(frame))
 


 
 %%
 %trying to calibrate the gain of the camera so that frames with different exposures can be
 %compared
figure(1)
exposure_list=min(exposure_lim):3:max(exposure_lim);
exposure_list=-11
%gain_list=linspace(min(gain_lim),max(gain_lim),50);
gain_list=linspace(0,100,500);
frame_struct=[];
cam.Brightness=0;
ii=1;

for m=1:size(exposure_list,2)
    cam.Exposure=exposure_list(m); 
    for n=1:size(gain_list,2)
        cam.Gain=gain_list(n);
        pause(0.01)
        frame_struct.setings(ii,:)=[cam.Exposure cam.Gain cam.Brightness];
        frame = double(snapshot(cam))/(2^8-1);
        frame=squeeze(frame(:,:,1));
        [N,edges] = histcounts(frame(:),linspace(0,1,1e5));
        centers=edges(1:end-1)+diff(edges)/2;
        frame_struct.hist(ii,:,:)=[N;centers];
        frame_struct.stat(ii,:)=[min(frame(:)),max(frame(:)),mean(frame(:)),std(frame(:))];
        %clf
        %subplot(2,1,1)
        %plot(centers,N.^0.1)
        %subplot(2,1,2)
        %imagesc(frame)
        %colormap('gray')
        fprintf('%2i,%2.1f,%2.2f\n',exposure_list(m),gain_list(n),frame_struct.stat(ii,2))
        ii=ii+1;
    end
end


%% lets try and guess the function that combines gain and exposure

mask=frame_struct.stat(:,2)<0.8;
exp_gain_comb=0.05178+(0.000806)*frame_struct.setings(mask,2);
max_norm=frame_struct.stat(mask,2)-frame_struct.stat(mask,3); %max-mean
poly=polyfit(frame_struct.setings(mask,2),max_norm,1);
%poly=[0.0000020371 0.0006152301 0.0543738513];
poly=[0.0008189398 0.0509854929 ];
%(frame_struct.setings(mask,1).^2);
figure(1)
clf;
set(gcf,'Color',[1 1 1]);
subplot(1,2,1)
plot(frame_struct.setings(mask,2),max_norm,'x')
hold on
plot(frame_struct.setings(mask,2),polyval(poly,frame_struct.setings(mask,2)))

subplot(1,2,2)
plot(frame_struct.setings(mask,2),(frame_struct.stat(mask,2)-frame_struct.stat(mask,3))./exp_gain_comb,'x')


 %%
 %trying to calibrate the gain of the camera so that frames with different exposures can be
 %compared
figure(1)
exposure_list=min(exposure_lim):3:max(exposure_lim);
exposure_list=-13:7;
%gain_list=linspace(min(gain_lim),max(gain_lim),50);
gain_list=linspace(0,100,50);
frame_struct=[];
cam.Brightness=0;
ii=1;

for m=1:size(exposure_list,2)
    cam.Exposure=exposure_list(m); 
    for n=1:size(gain_list,2)
        cam.Gain=gain_list(n);
        pause(0.01)
        frame_struct.setings(ii,:)=[cam.Exposure cam.Gain cam.Brightness];
        frame = double(snapshot(cam))/(2^8-1);
        frame=squeeze(frame(:,:,1));
        [N,edges] = histcounts(frame(:),linspace(0,1,1e5));
        centers=edges(1:end-1)+diff(edges)/2;
        frame_struct.hist(ii,:,:)=[N;centers];
        frame_struct.stat(ii,:)=[min(frame(:)),max(frame(:)),mean(frame(:)),std(frame(:))];
        %clf
        %subplot(2,1,1)
        %plot(centers,N.^0.1)
        %subplot(2,1,2)
        %imagesc(frame)
        %colormap('gray')
        fprintf('%2i,%2.1f,%2.2f\n',exposure_list(m),gain_list(n),frame_struct.stat(ii,2))
        ii=ii+1;
    end
end


%% lets try and guess the function that combines gain and exposure

mask=frame_struct.stat(:,2)<0.8;
%poly=[0.0000020371 0.0006152301 0.0543738513];
poly=[0.0008189398 0.0509854929 ];
exp_gain_comb=1800*polyval(poly,frame_struct.setings(mask,2)).*(1.8.^frame_struct.setings(mask,1));
max_norm=frame_struct.stat(mask,2)-frame_struct.stat(mask,3); %max-mean
%poly=polyfit(frame_struct.setings(mask,2),max_norm,1);
%(frame_struct.setings(mask,1).^2);
figure(1)
clf;
set(gcf,'Color',[1 1 1]);
subplot(1,2,1)
plot(frame_struct.setings(mask,2),max_norm,'x')
hold on
plot(frame_struct.setings(mask,2),exp_gain_comb,'rx')

subplot(1,2,2)
plot(frame_struct.setings(mask,2),(frame_struct.stat(mask,2)-frame_struct.stat(mask,3))./exp_gain_comb,'x')


%% couldnt find it in future the thing to do is to fit the peakval vs exposure curve for each exposure
%then use this to find a paramterization or even just table of the coefs




%%
idx=40;
clf
frame_struct.stat(idx,3)
plot(squeeze(frame_struct.hist(idx,2,:)),squeeze(frame_struct.hist(idx,1,:)).^0.25)

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
figure(1)
clf
plot(frame_struct.stat(:,2),'k')
hold on
plot(frame_struct.setings(:,2)/100,'r')


%% does brightness change their 

figure(1)

cam.Exposure=-13; 
snapshot(cam);
cam.Gain=50; 
snapshot(cam);
blist=-64:64;%-64:64%-64:0;
for bval=blist
pause(0.1)
cam.Brightness=bval;
snapshot(cam);
pause(0.1)

frame_struct.setings(ii,:)=[cam.Exposure cam.Gain cam.Brightness];
frame = double(snapshot(cam))/(2^8-1);
frame=squeeze(frame(:,:,1));
[N,edges] = histcounts(frame(:),linspace(0,1,1e5));
centers=edges(1:end-1)+diff(edges)/2;
plot(centers,N.^0.05)
end

%%

figure(1)

cam.Exposure=-2; 
snapshot(cam);
cam.Gain=50; 
pause(0.1)
snapshot(cam);
pause(0.1)
bval=-40;
cam.Brightness=bval;
cam.Brightness=bval;
cam.Sharpness=0;
pause(0.1)
snapshot(cam);
pause(0.1)


frame = double(snapshot(cam))/(2^8-1);
frame=squeeze(frame(:,:,1));
[N,edges] = histcounts(frame(:),linspace(0,1,1e5));
centers=edges(1:end-1)+diff(edges)/2;
plot(centers,N.^0.05)



%%
clear('cam');