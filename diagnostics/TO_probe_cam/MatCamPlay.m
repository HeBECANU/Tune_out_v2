%some scripts to play with http://www.elpcctv.com/full-hd-2mp-usb-camera-black-and-white-module-p-231.html
%eventualy want to build into a program to be used for beam profiling
%find the peak region, add lines to main plot
%subplot of this peak region with gaussian fit (line/2d)
%automatic exposure control and linearization to give a measure of power
%hilight saturated regions
%%
info=imaqhwinfo('winvideo');
formats=info.DeviceInfo.SupportedFormats;

vid = videoinput('winvideo', 1,'MJPG_1920x1080');
triggerconfig(vid,'manual')
addpath('Colormaps')

%%
src = getselectedsource(vid);
vid.FramesPerTrigger = 1;
preview(vid);
src.BacklightCompensation='off'
src.ExposureMode='manual';
src.Exposure=-1; %[-13 -1] %t=2^exposure
src.Contrast=0; %[0 64]
src.Sharpness=0;
src.Gain=10;
src.Gamma=100; %vout+A vin ^ 100/gamma [72 500]
src.Brightness=-30; %[-64 64]
src.Hue=0;
src.WhiteBalanceMode='manual';

%%

propinfo(src,'FrameRate')
exp=propinfo(src,'Exposure')
propinfo(src,'ExposureMode')


 %%
frame = getsnapshot(vid);
frame=squeeze(frame(:,:,1));
imagesc(flipud(frame'))
%colormap('veridis')
colormap(inferno())
caxis([0 intmax('uint8')])
pbaspect([1,1,1])
size(frame)
max(max(frame))
 
%% test frame rate
tic
times=zeros(200,1);
for n=1:size(times,1)
    src.Gain=mod(n,10)*10;
    frame = getsnapshot(vid);
    imagesc(frame)
    times(n)=toc;
    tic
end
mean(1./times)


%%

bdata=[];
b=-60:64;
for n=1:size(b,2)
    src.Brightness=b(n);
    %getsnapshot(vid);
    frame = getsnapshot(vid);
    frame=squeeze(frame(:,:,1));
    hist(single(frame(:)),100)
    xlim([0,255])
    %imagesc(frame)
    %colormap('gray')
    size(frame);
    bdata(n,1)=max(frame(:));
    bdata(n,2)=mean(frame(:));
    bdata(n,3)=std(single(frame(:)));
    pause(0.01)
end

 plot(b,bdata(:,1))
 set(gcf,'Color',[1 1 1]);
 hold on
  plot(b,bdata(:,2),'k')
  plot(b,bdata(:,3),'b')
 hold off

%%
src.Brightness=0;
frame = getsnapshot(vid);
frame=squeeze(frame(:,:,1));


 
 
 %%
expo=propinfo(src,'Exposure');
bright=propinfo(src,'Brightness');
gain=propinfo(src,'Gain');
figure(1)
bdata=[];
blim=bright.ConstraintValue;
elim=expo.ConstraintValue;
glim=gain.ConstraintValue;
b=min(blim):max(blim);
e=min(elim):max(elim);
src.Gain=0;
frame_struct=[];
frame = single(getsnapshot(vid))/(2^8-1);
frame=squeeze(frame(:,:,1));
triggerconfig(vid, 'manual');
ii=1;
for m=1:size(e,2)
    src.Exposure=e(m); 
    for n=1:size(b,2)
        tic
        src.Brightness=b(n);
        fprintf('%2i,%2i \n',b(n),e(m))
        frame_struct.setings(ii,:)=[src.Exposure src.Gain src.Brightness];
        
        
        frame =double(peekdata(vid,1));
        %frame = double(getsnapshot(vid))/(2^8-1);
        frame=squeeze(frame(:,:,1));
        [N,edges] = histcounts(frame(:),linspace(0,1,1e5));
        
        centers=edges(1:end-1)+diff(edges)/2;
        frame_struct.hist(ii,:)=[N,centers];
        frame_struct.stat(ii,:)=[min(frame(:)),max(frame(:)),mean(frame(:)),std(frame(:))];
        %clf
        %subplot(2,1,1)
        %plot(centers,N.^0.1)
        %subplot(2,1,2)
        %imagesc(frame)
        %colormap('gray')
        %pause(0.1)
        ii=ii+1;
    end
end

 set(gcf,'Color',[1 1 1]);
 hold on
  plot(b,bdata(:,2),'k')
  plot(b,bdata(:,3),'b')
 hold off


%%


expo=propinfo(src,'Exposure');
bright=propinfo(src,'Brightness');

figure(2)
b=-64;
blim=bright.ConstraintValue;
elim=expo.ConstraintValue;
target_maxval=0.5;
src.Exposure=-13; 
blur=5;
thesh=0.3;
area_thresh=10;

while true
    src.Brightness=b;
    frame_sd = double(getsnapshot(vid))./double(intmax('uint8')); %scaled and in double format
    frame_sd=double(squeeze(frame_sd(:,:,1)));
    frame_sd=flipud(frame_sd');
    frame_smooth=imgaussfilt(frame_sd,blur,'Padding','replicate');
    stat=[max(frame_sd(:)),mean(frame_sd(:)),std(single(frame_sd(:)))];
    frame_thresh=frame_smooth>thesh;
    regions=regionprops('table',frame_thresh);
    regions=sortrows(regions,1,'descend');
    areas=regions(:,:).Area;
    area_mask=areas>area_thresh;
    centers = regions(:,:).Centroid;
    
    
    subplot(1,4,1)
    imagesc(frame_sd)
    colormap(inferno())
    caxis([0 1])
    pbaspect([1,1,1])

    subplot(1,4,2)
    imagesc(frame_smooth)
    colormap(inferno())
    caxis([0 1])
    pbaspect([1,1,1])

    subplot(1,4,3)
    imagesc(frame_thresh)
    colormap(inferno())  
    if ~isequal(centers,[])
        hold on
        scatter(centers(area_mask,1),centers(area_mask,2),100,'r','+')
        hold off
    
    
    center=round(regions.Centroid(1,:));
    radius=round(regions.BoundingBox(1,3:4))*2;
    minpos=center-radius;
    minpos=max(1,minpos);
    maxpos=center+radius;
    maxpos=min(maxpos,fliplr(size(frame_sd)));

    subplot(1,4,4)
    imagesc(frame_sd(minpos(2):maxpos(2),minpos(1):maxpos(1)))

    else
    
        subplot(1,4,4)
        imagesc(zeros(1))
        caxis([0 1])
    
    end
    
    fprintf('%3i,%3i\n',stat(1),b)
    if stat(1)>target_maxval && b>blim(1)
        b=b-1;
    elseif stat(1)<target_maxval && b<blim(2)
        b=b+1;
    end
   
    pause(0.01)
end
fprintf('\n')



%%
delete(vid);