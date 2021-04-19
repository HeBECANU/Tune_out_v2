function out_st=parse_folder_name_into_wp_angles(path)

fname=split(fileparts(path),filesep);
fname=fname(end);
%% path reserved for future use
fname_lowercase=lower(fname);
fname_split=split(fname_lowercase,'_');
hwp_idx=strcmp(fname_split,'hwp');
if sum(hwp_idx)==0  
    hwp_ang=nan;
elseif sum(hwp_idx)>1
      error('multiple hwp matches')  
else
    hwp_ang=str2double(fname_split(find(hwp_idx,1)+1));
end


qwp_idx=strcmp(fname_split,'qwp');
if sum(qwp_idx)==0  
    qwp_ang=nan;
elseif sum(qwp_idx)>1
      error('multiple hwp matches')  
else
    qwp_ang=str2double(fname_split(find(qwp_idx,1)+1));
end

if isnan(hwp_ang) && isnan(qwp_ang)
    warning('no waveplate')
end

out_st=[];
out_st.qwp_ang=qwp_ang;
out_st.hwp_ang=hwp_ang;


end