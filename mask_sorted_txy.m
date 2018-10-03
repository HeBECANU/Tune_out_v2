function txymasked=mask_sorted_txy(txyin,txylim)
%ASSUMES TIME SORTED DATA! must pass  issorted(txyin(:,1))
%a fast masking function that gets used a lot in our data processing
%txylim=[[tmin,tmax];[xmin,xmax];[ymin,ymax]] (in seconds and meters)
%uses binary compare in order to give a speedup of approx 2 with LARGE 3e6+ count datasets

%2018-10-01
%   -speedup seems to be limited by the overheads in calling fast_sorted_mask
%   -perhaps compilation could overcome this
tvec=txyin(:,1);
mask=fast_sorted_mask(tvec,txylim(1,1),txylim(1,2));
txymasked=txyin(mask(1):mask(2),:);
mask_xy_subest=txymasked(:,2)>txylim(2,1) & txymasked(:,2)<txylim(2,2) &...
     txymasked(:,3)>txylim(3,1) & txymasked(:,3)<txylim(3,2);
txymasked=txymasked(mask_xy_subest,:);         
end
