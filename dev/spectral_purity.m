%filt_to = cell(1,4)'
filt_to=[0,725735894.6,40,725735894.6,1000;
        1,725736039.0,74,725735834.9,24;
        2,725735834.1,61,725735851.8,33;
        3,725735811.4,73,725735906.5,34;
    ];
filt_to_avg = zeros(4,3);
for jj = 1:size(filt_to,1)
    filt_to_avg(jj,:) = [filt_to(jj,1),(filt_to(jj,2)./filt_to(jj,3).^2+filt_to(jj,4)./filt_to(jj,5).^2)./(1./filt_to(jj,3).^2+1./filt_to(jj,5).^2),sqrt(var(filt_to(jj,[2,4]),1./filt_to(jj,[3,5]).^2,2))];
end
figure(451)
set(gcf,'color','w')
errorbar(filt_to_avg(:,1),filt_to_avg(:,2)-mean(filt_to_avg(:,2)),filt_to_avg(:,3),'kx')

%%
filt_to=[0,725735894.6,40;
        1,725736039.0,74;
        1,725735834.9,24;
        1,725735897.6,35;
        2,725735834.1,61;
        2,725735851.8,33;
        3,725735811.4,73;
        3,725735906.5,34;
    ];
figure(452)
set(gcf,'color','w')
errorbar(filt_to(:,1),filt_to(:,2)-mean(filt_to(:,2)),filt_to(:,3),'kx')
xlim([-0.5,3.5])