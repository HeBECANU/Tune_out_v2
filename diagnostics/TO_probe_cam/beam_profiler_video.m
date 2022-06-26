%%
clc; clear;
fclose('all')
objects = imaqfind;
delete(objects);


%%
addpath('Colormaps')
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
cam.Exposure=-10; %[-13 -1] %t=2^exposure
cam.Contrast=0; %[0 64]
cam.Sharpness=0;
cam.Gamma=100; %vout+A vin ^ 100/gamma [72 500]
cam.Brightness=-10; %[-64 64]
cam.Hue=0;
cam.WhiteBalanceMode='manual';
cam.Gain=73;

src=getselectedsource(vid);
frame_rate=src.FrameRate;


%webcam has no method for allowed values so i these will be hard coded for now
property=propinfo(cam,'Exposure');
exposure_lim=property.ConstraintValue;
brightness_lim=[-64 64];
exposure_lim=[-13 -2];
gamma_lim=[0 100];
gain_lim=[0 100];

%%
clear('imageData','time_stamp');
fprintf('starting to acquire video\n')
%triggerconfig(vid, 'Manual')
vid.FramesPerTrigger = inf;
start(vid)
%trigger(vid)
pause(60)
stop(vid)
numAvail = vid.FramesAvailable;
[imageData, time_stamp]  = getdata(vid,numAvail);
fprintf('video aquired\n')
date_string=datestr(datetime('now'),'yyyymmddTHHMMSS.FFF')

%%
%write to AVI slow as fuck
fprintf('saving video as avi\n')
v = VideoWriter(['./videos/',date_string,sprintf('_%.2fs_%.3fHz',[size(imageData,4)*mean(diff(time_stamp)),1/mean(diff(time_stamp))]),'.avi']);
v.FrameRate=1/mean(diff(time_stamp));
open(v)
writeVideo(v,imageData)
close(v)
fprintf('done video saved\n')

%%
fprintf('saving video\n')
save(['./videos/',date_string,'.mat'],'imageData','time_stamp','-v7.3')
fprintf('done video saved\n')


