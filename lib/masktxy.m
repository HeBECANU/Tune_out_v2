function txymasked=masktxy(txyin,txylim)
%simple masking function that gets used a lot in our data processing
%txylim=[[tmin,tmax];[xmin,xmax];[ymin,ymax]] (in seconds and meters)
mask=txyin(:,1)>txylim(1,1) & txyin(:,1)<txylim(1,2) &...
     txyin(:,2)>txylim(2,1) & txyin(:,2)<txylim(2,2) &...
     txyin(:,3)>txylim(3,1) & txyin(:,3)<txylim(3,2);
txymasked=txyin(mask,:);
            
end