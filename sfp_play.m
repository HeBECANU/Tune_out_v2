path='\\amplpc29\Users\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\2018-09-16_half_wp_angle_328\log_sfp_20180915T012012.txt';
fid = fopen(path,'r');
sfp_log_file_cells=textscan(fid,'%s','Delimiter','\n');
fclose(fid);


%%
