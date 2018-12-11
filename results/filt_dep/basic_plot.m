
filt_dep=[[3,362867189.6,48];[0,362867361.8,114];[0,362867618.6,60];[2,362867582.2,34];[2,362867550.1,26]];
mean_with_three=nanmean(filt_dep(filt_dep(:,1)==3,2));
figure(3)
errorbar(filt_dep(:,1),filt_dep(:,2)-mean_with_three,filt_dep(:,3),'xk')
xlabel('Number of filters')
ylabel('Change in TO (MHz)')
xlim([-0.5,3.5])

saveas(gcf,'.\results\filt_dep\plot.png')



%%
%last filter skew
%[cen filter,measured To,unc]
%aom correction not applied to any of this data
to_value=362867621*2;
filt_dep=[[to_value-111000,362867661.4*2,66*2];...  %quad Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20181201_filt_skew_neg111ghz\out\20181203T114556
    [to_value+110000,362867689.7*2,89*2];...  %quad Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20181202_filt_skew_pos110ghz\out\20181203T091828
    [to_value+50000,362867443.9*2,58*2];...  %quad Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20181203_filt_skew_pos50ghz_bad_setpt\out\20181203T135012
    [to_value+50000,725734744.4,331];... %quad  %Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20181203_filt_skew_pos50ghz\out\20181203T195730
    [to_value-50000,725736090.2,157];... %quad %Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20181202_filt_skew_neg50ghz\out\20181203T150450
    [to_value,725735356.2,187];... %quad  Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20181123_3_filt_align_dep_31um\out\20181203T203502
    [to_value,362867189.6*2,44*2];...  %Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20181120_filt_dep_3filt\out\20181121T220758
    [to_value-50000,725735786.3,115];... %quad Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20181204_filt_skew_neg50ghz\out\20181204T112057
    [to_value,725735093.4,102];...  %quad Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20181204_baseline_1\out\20181204T135302
    [to_value,725734668.0,179];... %quad Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20181120_filt_dep_3filt\out\20181204T153114
    [to_value,725734365.8,135]];
%inclue data from

figure(3)
errorbar((filt_dep(:,1)-to_value)*1e-3,filt_dep(:,2)-to_value,filt_dep(:,3),'xk')
xlabel(sprintf('Filter Center (GHZ) - %.1f (MHz) (blue)',to_value))
ylabel(sprintf('Measured To - %.1f (MHz) (blue) ',to_value))
xlim([-250,250])