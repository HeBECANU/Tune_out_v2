anal_opts.dir='Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20181818_power_stab_2';
csv_files=dir(fullfile(anal_opts.dir,'*.csv'));
file_name=fullfile(anal_opts.dir,csv_files.name);
power_dat=importpowermeter(file_name);
power_val_times=[power_dat.time,power_dat.data];

power_thresh=8e-3;
mask=power_val_times(:,2)>power_thresh;
figure(1)
set(gcf,'color','w')
plot(power_val_times(mask,1),power_val_times(mask,2)*1e3)
xlabel('Time(s)')
ylabel('Power (mW)')


