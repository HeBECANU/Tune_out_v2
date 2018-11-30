
filt_dep=[[3,362867189.6,48];[0,362867361.8,114];[0,362867618.6,60];[2,362867582.2,34];[2,362867550.1,26]];
mean_with_three=nanmean(filt_dep(filt_dep(:,1)==3,2));
figure(3)
errorbar(filt_dep(:,1),filt_dep(:,2)-mean_with_three,filt_dep(:,3),'xk')
xlabel('Number of filters')
ylabel('Change in TO (MHz)')
xlim([-0.5,3.5])

saveas(gcf,'.\results\filt_dep\plot.png')