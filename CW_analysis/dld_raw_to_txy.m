%%% Use this to convert raw files from DLD to txy format
%%% which is used by C++ program such as
%%% g2_calc_norm_across_files_spatial_dy.exe to compute g2
%%% there is a 1:1 correspondance between input files and 
%%% output files with '_txy_forc'


function[] = dld_raw_to_txy(filename_raw,start_file,end_file)
%tic

for i=start_file:end_file
    
    file_no = num2str(i);
    filename_read = [filename_raw,file_no];
    
    [hits_sorted] = dld_read_5channels_reconst_multi_imp(filename_read,1,0,1,0);
    
    filename_write = [filename_raw,'_txy_forc',file_no,'.txt'];
    
    dlmwrite(filename_write,hits_sorted,'precision',8);
   
    %toc
end